"""Synthetic people and their timed daily movement."""

from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace
from itertools import groupby, pairwise
from math import atan2
from numbers import Integral, Real

import numpy as np
from shapely import line_locate_point, points
from shapely.geometry import LineString, Point, box

from ..geometry import Bounds
from .object import Object
from .roadnetwork import RoadNetwork

DAY_SECONDS = 86_400
STATES = frozenset({"indoor", "outdoor", "wait", "walk", "cycle", "car", "transit"})


def _walking_stages(paths):
    stages = []
    for path in paths:
        walking = path.state in {"walk", "wait"}
        if stages:
            previous = stages[-1]
            if (
                not walking
                or path.journey_id != previous.journey_id
                or path.start != previous.end
                or not np.array_equal(path.samples[0, 1:], previous.samples[-1, 1:])
            ):
                yield stages
                stages = []
        if walking:
            stages.append(path)
    if stages:
        yield stages


def _edge_runs(paths):
    """Join walking stages and retain an intervening wait on a continuing edge."""
    for stages in _walking_stages(paths):
        edges = []
        for path in stages:
            if path.state == "wait":
                edges.extend([None] * (len(path.samples) - 1))
            else:
                edges.extend(path.edge_ids or ("",) * (len(path.samples) - 1))
        runs = [(edge, len(list(steps))) for edge, steps in groupby(edges)]
        for i in range(1, len(runs) - 1):
            edge, length = runs[i]
            before, after = runs[i - 1][0], runs[i + 1][0]
            if edge is None and before and before == after:
                runs[i] = before, length
        if not any(edge for edge, _ in runs):
            continue
        samples = stages[0].samples
        if len(stages) > 1:
            samples = np.concatenate([samples, *(p.samples[1:] for p in stages[1:])])
        start = 0
        for edge, pieces in groupby(runs, key=lambda run: run[0]):
            stop = start + sum(length for _, length in pieces)
            if edge:
                yield edge, samples[start : stop + 1]
            start = stop


def _edge_line(xy, *, inferred=False):
    """Orient open edges lexicographically and preserve the network order of loops."""
    start, end = tuple(xy[0]), tuple(xy[-1])
    if start > end or (
        inferred
        and start == end
        and tuple(map(tuple, xy[::-1])) < tuple(map(tuple, xy))
    ):
        xy = xy[::-1]
    return LineString(xy[np.r_[True, np.any(xy[1:] != xy[:-1], axis=1)]])


def _edge_traversals(samples, line=None):
    """Split a timed run whenever the pedestrian reverses along the edge."""
    xy = samples[:, 1:3]
    steps = np.diff(xy, axis=0)
    lengths = np.linalg.norm(steps, axis=1)
    moving = np.flatnonzero(lengths)
    if not len(moving):
        return
    unit = steps[moving] / lengths[moving, None]
    signs = np.zeros(len(moving))
    if line is not None:
        progress = np.diff(line_locate_point(line, points(xy)))[moving]
        if line.is_closed:
            length = line.length
            progress[progress > length / 2] -= length
            progress[progress < -length / 2] += length
        signs = np.sign(np.where(abs(progress) > 1e-8, progress, 0.0))
    start = 0
    direction = signs[0]
    for k, i in enumerate(moving[1:], 1):
        before, after = signs[k - 1 : k + 1]
        if before and after:
            reversing = before != after
        else:
            # Compare the sine of the turn angle independently of step length.
            u, v = unit[k - 1 : k + 1]
            reversing = u @ v < 0 and abs(u[0] * v[1] - u[1] * v[0]) <= 1e-6
        if reversing:
            yield samples[start : i + 1], direction
            start = i
            direction = after
        elif not direction:
            direction = after
    yield samples[start:], direction


def _screenline_endpoints(start, end):
    a, b = np.asarray(start, dtype=float), np.asarray(end, dtype=float)
    if a.shape != (2,) or b.shape != (2,):
        raise ValueError("Screenline requires distinct finite XY endpoints")
    direction = b - a
    length_squared = float(direction @ direction)
    if not np.isfinite(length_squared) or not length_squared:
        raise ValueError("Screenline requires distinct finite XY endpoints")
    return a, direction, length_squared


def _screenline_contact(samples, signed, a, direction, length_squared):
    """Return the first contact and its incoming segment between opposite sides."""
    if len(samples) == 2:
        fraction = signed[0] / (signed[0] - signed[1])
        point = samples[0] + fraction * (samples[1] - samples[0])
        along = (point[1:3] - a) @ direction / length_squared
        return (point, 0) if 0 < along < 1 else (None, None)
    along = (samples[1:-1, 1:3] - a) @ direction / length_squared
    for i, pos in enumerate(along, 1):
        if 0 < pos < 1:
            return samples[i], i - 1
        if i < len(along):
            following = along[i]
            if (pos <= 0 < following) or (pos >= 1 > following):
                endpoint = 0 if pos <= 0 else 1
                fraction = (endpoint - pos) / (following - pos)
                return samples[i] + fraction * (samples[i + 1] - samples[i]), i
    return None, None


def _numeric_array(values):
    array = np.asarray(values, dtype=object)
    if any(
        isinstance(value, bool) or not isinstance(value, Real) for value in array.flat
    ):
        raise ValueError("coordinates and samples must be numeric")
    return array.astype(float)


@dataclass(frozen=True)
class MovementState:
    """Interpolated location, heading in radians clockwise from north, and speed in m/s."""

    x: float
    y: float
    z: float
    heading: float
    speed: float
    state: str
    purpose: str
    distance: float


@dataclass(eq=False, frozen=True)
class TimedPath:
    """An activity or travel stage over [start, end) in local clock seconds.

    Samples have columns (second, easting, northing, height). Edge IDs identify
    segments between successive samples. An empty edge ID denotes an off-network
    connector. Stationary intervals repeat the same position at their endpoints.
    """

    state: str
    purpose: str
    samples: np.ndarray
    place_id: str = ""
    journey_id: str = ""
    edge_ids: tuple[str, ...] = ()

    def __post_init__(self):
        samples = np.array(self.samples, dtype=np.float64, order="C", copy=True)
        if (
            samples.ndim != 2
            or samples.shape[1] != 4
            or len(samples) < 2
            or not np.isfinite(samples).all()
            or np.any(np.diff(samples[:, 0]) <= 0)
            or samples[0, 0] < 0
            or samples[-1, 0] > DAY_SECONDS
        ):
            raise ValueError(
                "TimedPath requires strictly timed finite (second, x, y, z) samples within one day"
            )
        if (
            self.state not in STATES
            or not isinstance(self.purpose, str)
            or not self.purpose
        ):
            raise ValueError("TimedPath requires a supported state and a purpose")
        for name in ("place_id", "journey_id"):
            if not isinstance(getattr(self, name), str):
                raise ValueError(f"TimedPath {name} must be a string")
        edge_ids = tuple(self.edge_ids)
        if edge_ids and (
            len(edge_ids) != len(samples) - 1
            or not all(isinstance(e, str) for e in edge_ids)
        ):
            raise ValueError("TimedPath edge IDs must be strings aligned to segments")
        if self.state in {"indoor", "outdoor", "wait"} and np.any(
            samples[:, 1:] != samples[0, 1:]
        ):
            raise ValueError("An activity or wait must have a stationary location")
        samples.setflags(write=False)
        object.__setattr__(self, "samples", samples)
        object.__setattr__(self, "edge_ids", edge_ids)

    def __deepcopy__(self, memo):
        # The constructor preserves the read-only sample array contract.
        result = replace(self)
        memo[id(self)] = result
        return result

    def __eq__(self, other):
        return (
            isinstance(other, TimedPath)
            and (
                self.state,
                self.purpose,
                self.place_id,
                self.journey_id,
                self.edge_ids,
            )
            == (
                other.state,
                other.purpose,
                other.place_id,
                other.journey_id,
                other.edge_ids,
            )
            and np.array_equal(self.samples, other.samples)
        )

    @property
    def start(self) -> float:
        """Inclusive start in seconds from local midnight."""
        return float(self.samples[0, 0])

    @property
    def end(self) -> float:
        """Exclusive end in seconds from local midnight."""
        return float(self.samples[-1, 0])

    @property
    def length(self) -> float:
        """Path length in metres."""
        return float(np.linalg.norm(np.diff(self.samples[:, 1:], axis=0), axis=1).sum())

    def position_at(self, second: float) -> MovementState | None:
        """Interpolate on the authoritative segment, preserving corners and waits."""
        if not self.start <= second < self.end:
            return None
        idx = int(np.searchsorted(self.samples[:, 0], second, side="right")) - 1
        before, after = self.samples[idx : idx + 2]
        delta = after[1:] - before[1:]
        fraction = (second - before[0]) / (after[0] - before[0])
        distance = float(np.linalg.norm(delta))
        travelled = np.linalg.norm(
            np.diff(self.samples[: idx + 1, 1:], axis=0), axis=1
        ).sum()
        xyz = before[1:] + fraction * delta
        return MovementState(
            *map(float, xyz),
            atan2(delta[0], delta[1]),
            distance / (after[0] - before[0]),
            self.state,
            self.purpose,
            float(travelled + fraction * distance),
        )


@dataclass(repr=False)
class Person(Object):
    """Synthetic individual with stable scenario ID, demographics and daily stages."""

    household_id: str = ""
    age_band: str = "unknown"
    role: str = "unknown"
    home_id: str = ""
    regular_destination_id: str = ""
    appearance_key: int = 0
    paths: tuple[TimedPath, ...] = field(default_factory=tuple)

    def position_at(self, second: float) -> MovementState | None:
        """Return the person's state at local clock seconds or None outside their day."""
        index = bisect_right(self.paths, second, key=lambda path: path.start) - 1
        return self.paths[index].position_at(second) if index >= 0 else None

    @property
    def home_location(self) -> tuple[float, float, float] | None:
        """Return the inferred home anchor, including homes outside the view."""
        for path in self.paths:
            if path.purpose == "home" and path.state == "indoor":
                return tuple(map(float, path.samples[0, 1:]))
        return None

    def _validate_day(self):
        if (
            isinstance(self.appearance_key, bool)
            or not isinstance(self.appearance_key, Integral)
            or not 0 <= self.appearance_key < 2**32
        ):
            raise ValueError("Person appearance key must fit uint32")
        if any(type(path) is not TimedPath for path in self.paths):
            raise ValueError("Person daily paths must be TimedPath values")
        for before, after in pairwise(self.paths):
            if before.end != after.start or not np.array_equal(
                before.samples[-1, 1:], after.samples[0, 1:]
            ):
                raise ValueError(f"Person {self.id} has discontinuous daily paths")


@dataclass(frozen=True)
class PedestrianPassage:
    """One directional passage at a simulated local time."""

    person_id: str
    second: float
    direction: int
    edge_id: str


@dataclass(repr=False)
class PeopleSimulation(Object):
    """A local view of a synthetic day, with absolute coordinates in a projected metre CRS."""

    web_geojson = True

    def _summary_items(self):
        return [
            ("people", len(self.people)),
            ("reference_year", self.attributes.get("reference_year")),
            ("day_type", self.attributes.get("day_type")),
        ]

    def calculate_bounds(self, lod=None):
        """Return the simulation delivery extent in the projected CRS."""
        return self.view_bounds

    @property
    def people(self) -> list[Person]:
        """Return people associated with the view, including visitors with outside homes."""
        return self.get_children(Person)

    @property
    def view_bounds(self) -> Bounds:
        """Delivery area in the simulation CRS."""
        return Bounds(*self.attributes.get("view_bounds", ()))

    def _coordinate_reference(self):
        from pyproj import CRS
        from pyproj.exceptions import CRSError

        try:
            crs = CRS.from_user_input(self.transform.srs)
        except CRSError as error:
            raise ValueError(
                "PeopleSimulation requires a projected CRS in metres"
            ) from error
        if not crs.is_projected or any(
            axis.unit_conversion_factor != 1.0 for axis in crs.axis_info[:2]
        ):
            raise ValueError("PeopleSimulation requires a projected CRS in metres")
        if not np.array_equal(self.transform.affine, np.eye(4)):
            raise ValueError(
                "PeopleSimulation requires absolute path coordinates and an identity transform"
            )
        identity = np.eye(4)
        for person in self.people:
            frame = person.transform
            if not np.array_equal(frame.affine, identity) or (
                frame.srs and frame.srs != self.transform.srs and crs != frame.srs
            ):
                raise ValueError(
                    "Person paths must be in the simulation coordinate frame"
                )
        return crs

    def _vehicle_tracks(self):
        tracks = self.attributes.get("vehicles", ())
        if not isinstance(tracks, (tuple, list)):
            raise ValueError(
                "PeopleSimulation vehicles must be a list of scheduled tracks"
            )
        ids = set()
        for track in tracks:
            if not isinstance(track, dict) or any(
                not isinstance(track.get(field), str)
                for field in ("id", "kind", "line")
            ):
                raise ValueError(
                    "Scheduled vehicles require string id, kind and line fields"
                )
            vehicle_id = track["id"]
            if not vehicle_id.strip() or not track["kind"].strip() or vehicle_id in ids:
                raise ValueError(
                    "Scheduled vehicles require a unique nonempty id and a kind"
                )
            ids.add(vehicle_id)
            try:
                samples = TimedPath(
                    "transit", "service", _numeric_array(track.get("samples"))
                ).samples
            except (TypeError, ValueError, OverflowError) as error:
                raise ValueError(
                    f"Scheduled vehicle {vehicle_id!r}: {error}"
                ) from error
            yield track, samples

    def _screenline(self):
        line = self.attributes.get("screenline")
        if line is None:
            return None
        if not isinstance(line, dict) or not isinstance(line.get("name"), str):
            raise ValueError("Screenline requires a name")
        bin_seconds = line.get("bin_seconds")
        if type(bin_seconds) is not int or bin_seconds <= 0:
            raise ValueError("Screenline bin_seconds must be a positive integer")
        start, end = (_numeric_array(line.get(key)) for key in ("start", "end"))
        _screenline_endpoints(start, end)
        return line

    def _buildings(self):
        buildings = self.attributes.get("buildings", ())
        if not isinstance(buildings, (list, tuple)):
            raise ValueError("PeopleSimulation buildings must be a list")
        for building in buildings:
            if not isinstance(building, dict) or not isinstance(
                building.get("id"), str
            ):
                raise ValueError("Playback buildings require a string id")
            try:
                ring = _numeric_array(building.get("ring"))
            except (TypeError, ValueError, OverflowError) as error:
                raise ValueError(
                    "Playback buildings require a finite XY ring"
                ) from error
            if (
                ring.ndim != 2
                or ring.shape[1] != 2
                or not np.isfinite(ring).all()
                or len(np.unique(ring, axis=0)) < 3
            ):
                raise ValueError("Playback buildings require a finite XY ring")
            yield building, ring

    def _validate_playback(self):
        self._screenline()
        for _ in self._buildings():
            pass
        for _ in self._vehicle_tracks():
            pass

    def household_members(self, household_id: str) -> tuple[str, ...]:
        """Resolve membership even when some members or the household home are outside the view."""
        return tuple(self.attributes.get("households", {}).get(household_id, ()))

    def positions_at(self, second: float, *, outdoor_only: bool = False):
        """Iterate people and their states within the half-open delivery bounds."""
        bounds = self.view_bounds
        for person in self.people:
            state = person.position_at(second)
            if state is None or (
                outdoor_only and state.state not in {"walk", "wait", "outdoor"}
            ):
                continue
            if (
                bounds.xmin <= state.x < bounds.xmax
                and bounds.ymin <= state.y < bounds.ymax
            ):
                yield person, state

    def passages(
        self, start, end, *, interval=(0, DAY_SECONDS), edge_ids=None, height_range=None
    ):
        """Count oriented crossings of the open screenline from start to end.

        A transition from its left to right is +1. A walk along the line is
        timed at its first entry into the screenline segment. The clock
        interval is half-open. A touch and return contributes zero. Filter
        network edge IDs or heights when the line intersects a different level.
        """
        a, direction, length_squared = _screenline_endpoints(start, end)
        lower, upper = interval
        if not 0 <= lower <= upper <= DAY_SECONDS:
            raise ValueError(
                "Passage interval must be ordered within the simulation day"
            )
        allowed = None if edge_ids is None else set(edge_ids)
        events = []
        for person in self.people:
            for walking, stages in groupby(
                person.paths, lambda p: p.state in {"walk", "wait"}
            ):
                if not walking:
                    continue
                paths = list(stages)
                if paths[-1].end < lower or paths[0].start >= upper:
                    continue
                samples = np.concatenate(
                    [paths[0].samples, *(p.samples[1:] for p in paths[1:])]
                )
                edges = [
                    edge
                    for p in paths
                    for edge in (p.edge_ids or ("",) * (len(p.samples) - 1))
                ]
                relative = samples[:, 1:3] - a
                signed = direction[0] * relative[:, 1] - direction[1] * relative[:, 0]
                indices = np.flatnonzero(abs(signed) > 1e-8)
                for i, j in pairwise(indices):
                    if (signed[i] > 0) == (signed[j] > 0):
                        continue
                    point, segment = _screenline_contact(
                        samples[i : j + 1],
                        signed[i : j + 1],
                        a,
                        direction,
                        length_squared,
                    )
                    if point is None:
                        continue
                    edge = edges[i + segment]
                    if allowed is not None and edge not in allowed:
                        continue
                    if lower <= point[0] < upper and (
                        height_range is None
                        or height_range[0] <= point[3] < height_range[1]
                    ):
                        events.append(
                            PedestrianPassage(
                                person.id,
                                float(point[0]),
                                1 if signed[i] > 0 else -1,
                                edge,
                            )
                        )
        return sorted(events, key=lambda event: (event.second, event.person_id))

    def _road_edge_lines(self, requested):
        indexed, named = {}, {}
        for roads in self.get_children(RoadNetwork):
            count = len(roads.edges)
            columns = {
                name: values
                for name in ("id", "osm_way_id", "segment_index")
                if isinstance(
                    values := roads.attributes.get(name), (list, tuple, np.ndarray)
                )
                and (not isinstance(values, np.ndarray) or values.ndim == 1)
                and len(values) == count
            }
            ids = columns.get("id")
            ways = columns.get("osm_way_id")
            segment_indices = columns.get("segment_index", range(count))
            geometry = roads.linestrings
            for i, nodes in enumerate(roads.edges):
                edge_ids = []
                if ids is not None and ids[i] is not None:
                    edge_ids.append(str(ids[i]))
                if ways is not None and ways[i] is not None:
                    edge_ids.append(f"{ways[i]}:{segment_indices[i]}")
                keys = [
                    key for key in edge_ids if key in requested and key not in named
                ]
                idx = str(i)
                by_index = idx in requested and idx not in indexed
                if not keys and not by_index:
                    continue
                vertices = (
                    geometry[i].vertices
                    if len(geometry) == count
                    else roads.vertices[nodes]
                )
                xy = vertices[:, :2]
                if len(xy) >= 2 and np.any(xy != xy[0]):
                    line = _edge_line(xy)
                    named.update(dict.fromkeys(keys, line))
                    if by_index:
                        indexed[idx] = line
        return indexed | named

    def segment_passages(self, bin_seconds: int = 300):
        """Return directional per-edge counts and endpoints keyed by clock bin.

        A traversal is timed at its wall-clock midpoint, including waits on
        the continuing edge. Only midpoints in the delivery area contribute.
        Direction is positive toward the larger lexicographic XY endpoint,
        or along the network geometry on a loop. Without matching network
        geometry, the longest observed traversal supplies the edge reference.
        """
        if not isinstance(bin_seconds, int) or bin_seconds <= 0:
            raise ValueError("bin_seconds must be a positive integer")
        bounds = self.view_bounds
        counts = Counter()
        coordinates = {}
        runs = defaultdict(list)
        for person in self.people:
            for edge, samples in _edge_runs(person.paths):
                runs[edge].append(samples)
        road_lines = self._road_edge_lines(runs.keys())
        for edge, samples_list in runs.items():
            line = road_lines.get(edge)
            if line is None:
                candidates = (
                    _edge_line(traversal[:, 1:3], inferred=True)
                    for samples in samples_list
                    for traversal, _ in _edge_traversals(samples)
                )
                line = max(
                    candidates, key=lambda item: (item.length, item.wkb), default=None
                )
            if line is None:
                continue
            endpoints = line.coords[0], line.coords[-1]
            for samples in samples_list:
                for traversal, direction in _edge_traversals(samples, line):
                    a, b = traversal[0], traversal[-1]
                    second = (a[0] + b[0]) / 2
                    x = np.interp(second, traversal[:, 0], traversal[:, 1])
                    y = np.interp(second, traversal[:, 0], traversal[:, 2])
                    if not (
                        bounds.xmin <= x < bounds.xmax
                        and bounds.ymin <= y < bounds.ymax
                    ):
                        continue
                    if not direction:
                        direction = 1 if tuple(a[1:3]) < tuple(b[1:3]) else -1
                    direction = int(direction)
                    counts[int(second // bin_seconds), edge, direction] += 1
                    coordinates[edge, direction] = (
                        endpoints if direction > 0 else endpoints[::-1]
                    )
        return counts, coordinates

    def screenline_profile(self):
        """Return the configured demonstration screenline and its binned counts."""
        line = self._screenline()
        if line is None:
            return None
        counts = Counter(
            int(event.second // line["bin_seconds"])
            for event in self.passages(line["start"], line["end"])
        )
        return {**line, "counts": {str(key): count for key, count in counts.items()}}

    def plot(self, second: float = 28800, ax=None, *, animate=False, show=True):
        """Plot pedestrians at a local time or play the next three minutes at 10× speed."""
        from matplotlib.collections import PolyCollection

        from ...plotting.style import (
            DTCC_COLORS,
            apply_dtcc_style,
            get_axes,
            plot_line_segments,
            show_plot,
        )

        if ax is None:
            ax = get_axes()
            ax.ticklabel_format(style="plain", useOffset=False)
        roads = next(iter(self.get_children(RoadNetwork)), None)
        if roads is not None:
            plot_line_segments(
                ax, roads._plot_segments(), color="#4f5b62", linewidth=0.8
            )
        buildings = self.attributes.get("buildings", ())
        if buildings:
            ax.add_collection(
                PolyCollection(
                    [building["ring"] for building in buildings],
                    facecolor=DTCC_COLORS["dark_panel"],
                    edgecolor=DTCC_COLORS["muted"],
                    linewidth=0.5,
                )
            )
        apply_dtcc_style(
            ax, equal_aspect=True, xlabel="Easting (m)", ylabel="Northing (m)"
        )
        markers = ax.scatter(
            [],
            [],
            s=34,
            color=DTCC_COLORS["yellow"],
            edgecolors=DTCC_COLORS["dark"],
            linewidths=0.5,
            zorder=5,
        )
        line = self.attributes.get("screenline")
        if line is not None and not animate:
            ax.plot(
                *zip(line["start"], line["end"], strict=True),
                color="#ac672b",
                label=line["name"],
            )
            ax.legend()
        bounds = self.view_bounds
        ax.set(xlim=(bounds.xmin, bounds.xmax), ylim=(bounds.ymin, bounds.ymax))

        def draw(time):
            positions = [
                (p.x, p.y) for _, p in self.positions_at(time, outdoor_only=True)
            ]
            markers.set_offsets(np.asarray(positions).reshape(-1, 2))
            hour, remainder = divmod(int(time), 3600)
            minute, second_of_minute = divmod(remainder, 60)
            ax.set_title(
                f"Synthetic pedestrians · {hour:02}:{minute:02}:{second_of_minute:02}"
                f" · {len(positions)} visible"
            )

        draw(second)
        if animate:
            from matplotlib.animation import FuncAnimation

            ax._dtcc_people_animation = FuncAnimation(
                ax.figure,
                draw,
                frames=np.arange(second, min(second + 180, DAY_SECONDS)),
                interval=100,
                cache_frame_data=False,
            )
        show_plot(show)
        return ax

    def to_geojson(self, crs: str = "EPSG:4326") -> dict:
        """Derive timed browser trajectories in WGS84 from the native samples."""
        from pyproj import CRS, Transformer

        bounds = self.view_bounds
        view = box(*bounds.tuple)
        if CRS.from_user_input(crs) != CRS.from_epsg(4326):
            raise ValueError("People playback GeoJSON requires EPSG:4326")
        crs_transform = Transformer.from_crs(
            self._coordinate_reference(), crs, always_xy=True
        )
        features = []
        for person in self.people:
            age = person.attributes.get("synthetic_age")
            figure = {"child": True} if age is not None and age < 13 else {}
            for path in person.paths:
                if path.state not in {"walk", "wait", "indoor", "outdoor", "transit"}:
                    continue
                samples = path.samples
                coords = samples[:, 1:3]
                shape = (
                    Point(coords[0])
                    if np.all(coords == coords[0])
                    else LineString(coords)
                )
                if not view.intersects(shape):
                    continue
                x, y = crs_transform.transform(samples[:, 1], samples[:, 2])
                distance = np.r_[
                    0.0,
                    np.cumsum(np.linalg.norm(np.diff(samples[:, 1:], axis=0), axis=1)),
                ]
                features.append(
                    {
                        "type": "Feature",
                        "geometry": {
                            "type": "LineString",
                            "coordinates": np.column_stack(
                                (x, y, samples[:, 3])
                            ).tolist(),
                        },
                        "properties": {
                            "person": person.id,
                            "appearance": int(person.appearance_key),
                            "state": path.state,
                            "purpose": path.purpose,
                            "seconds": samples[:, 0].tolist(),
                            "metres": distance.tolist(),
                            **({"place": path.place_id} if path.place_id else {}),
                            **figure,
                        },
                    }
                )
        screenline = self.screenline_profile()
        if screenline is not None:
            for endpoint in ("start", "end"):
                screenline[endpoint] = list(
                    crs_transform.transform(*screenline[endpoint])
                )
        return {
            "type": "FeatureCollection",
            "features": features,
            "people": {
                "bounds_wgs84": list(crs_transform.transform_bounds(*bounds.tuple)),
                "screenline": screenline,
                "buildings": [
                    {
                        "id": building["id"],
                        "ring": np.column_stack(
                            crs_transform.transform(*ring.T)
                        ).tolist(),
                        **{
                            key: building[key]
                            for key in ("height", "height_source", "kind")
                            if key in building
                        },
                    }
                    for building, ring in self._buildings()
                ],
                "vehicles": [
                    {
                        "id": vehicle["id"],
                        "kind": vehicle["kind"],
                        "line": vehicle["line"],
                        "seconds": samples[:, 0].tolist(),
                        "coordinates": np.column_stack(
                            (
                                *crs_transform.transform(*samples[:, 1:3].T),
                                samples[:, 3],
                            )
                        ).tolist(),
                    }
                    for vehicle, samples in self._vehicle_tracks()
                ],
                **{
                    key: self.attributes[key]
                    for key in (
                        "reference_year",
                        "day_type",
                        "time_zone",
                        "showcase_start",
                        "daylight",
                    )
                    if key in self.attributes
                },
            },
        }
