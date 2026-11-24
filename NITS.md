# Nits

Nit findings from merged pull requests that were accepted but left unfixed.
The pull request that fixes one deletes its line. Keep entries sorted by path.

Entry format; the symbol and PR number are optional:

    - `path/to/file.py` (`symbol`): what is wrong; why it was left open. (#PR)

- `dtcc_core/datasets/publish.py` (`DatasetUploadClient.upload_package`): a token with a character outside latin-1 makes `http.client` raise `UnicodeEncodeError`, which escapes unwrapped instead of as a `DatasetUploadError`; its message holds only the position, and such a token cannot authenticate anyway, so it was left open.
- `dtcc_core/datasets/transport/vasttrafik.py` (`_fetch_token`, `_fetch_positions`): a failing live run under pytest's default or `-l` traceback prints these frames' credential arguments; the live workflow runs with `--tb=short`, which does not, so it was left open. (#139)
- `dtcc_core/model/exchange.py` (`_encode_geometry`): A `VolumeGrid` with a masked `z_levels` entry fails: a masked end level at save, a masked inner level at load. Only values assigned after construction can be masked. (#215)
- `dtcc_core/model/geometry/grid.py` (`VolumeGrid._validate_z_levels`): A mixed int and float list can round in `np.asarray` before the precision check, so it passes validation and the stored levels lose the difference. This needs integers beyond 2**53. (#215)
- `dtcc_core/model/values/dataset.py` (`DatasetValue.__eq__`): `numpy_array == DatasetValue(...)` returns NumPy's elementwise array, not one bool; `__array_ufunc__ = None` would fix it but also disables NumPy arithmetic on DatasetValue.
- `dtcc_core/plotting/face_map.py` (`_layers`): The sampling note can overlap the scale bar, mention absent context and omit that field faces were sampled. This affects annotation clarity. (#197)
- `dtcc_core/plotting/face_map.py` (`_limits`): An invalid range with one explicit limit does not identify the automatic bound in the error message. The invalid range is rejected. (#197)

# Deferred HTTP transport findings

These accepted limitations apply to the shared HTTP clients. Resolve them with
focused regression tests when extending the affected behavior.

| Finding | Current behavior and reason for deferral | Follow-up |
| --- | --- | --- |
| Malformed request headers | Backend validation can raise different exception types. Valid requests use the shared transport error contract; invalid-input normalization is deferred. | Validate malformed headers consistently before a network attempt. |
| Raw-body Content-Type | Requests and aiohttp choose different defaults when callers omit Content-Type. Callers can provide an explicit header. | Choose and test a shared default before a provider relies on implicit content typing. |
| Invalid CA paths | TLS setup can raise a configuration exception directly. This preserves the underlying configuration failure. | Define a shared configuration-error contract before normalizing these exceptions. |
| Duplicated sync and async loops | Shared policy owns budgets, retries and redirects, but the backend I/O loops remain separate. Refactoring them would expand the current change. | Keep conformance tests shared; consolidate only where it simplifies both implementations. |
| CA environment variables | Native and Requests clients honor REQUESTS_CA_BUNDLE and CURL_CA_BUNDLE even with trust_env=False; aiohttp uses the documented download SSL context. Explicit ca_bundle provides consistent trust configuration. | Decide on a shared environment policy before changing the documented behavior. |

Cache entries with a `Vary` response header are deliberately not stored. Before
adopting the transport in a provider that returns `Vary: Accept-Encoding`, decide
whether to support variant-aware caching; the current restriction is documented
in [HTTP transport](docs/io/http-transport.md).

The dependency lower-bound checker imports `packaging` through the locked
environment's transitive dependencies. If the checker must run in a minimal
environment, declare that dependency explicitly or remove the import.
