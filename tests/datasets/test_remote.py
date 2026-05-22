"""Tests for RemoteDatasetDescriptor and remote service registration."""

import pytest
from unittest.mock import patch, MagicMock

from dtcc_core.datasets.remote import RemoteDatasetDescriptor


class TestRemoteDatasetDescriptor:
    def test_construction(self):
        desc = RemoteDatasetDescriptor(
            name="test_dataset",
            description="A test remote dataset",
            args_schema={"type": "object", "properties": {"bounds": {"type": "array"}}},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="dtcc-sim",
            timeout_hint=600,
        )
        assert desc.name == "test_dataset"
        assert desc.description == "A test remote dataset"
        assert desc.source_service == "dtcc-sim"
        assert desc.data_category == "simulation"
        assert desc.timeout_hint == 600
        assert desc.supported_formats == ["xdmf"]
        assert desc.base_url == "http://localhost:8001"

    def test_base_url_trailing_slash_stripped(self):
        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001/",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="test",
        )
        assert desc.base_url == "http://localhost:8001"
        assert desc.data_category == "remote"

    def test_explicit_data_category_is_preserved(self):
        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="dtcc-sim",
            data_category="derived",
        )
        assert desc.data_category == "derived"

    def test_dtcc_sim_data_category_detection_is_case_insensitive(self):
        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="DTCC-SIM",
        )
        assert desc.data_category == "simulation"

    def test_does_not_auto_register(self):
        """register=False prevents __init_subclass__ auto-registration."""
        from dtcc_core.datasets import list as list_datasets

        before = set(list_datasets().keys())
        _ = RemoteDatasetDescriptor(
            name="should_not_register",
            description="test",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="test",
        )
        after = set(list_datasets().keys())
        assert "should_not_register" not in after - before

    def test_show_options_returns_stored_schema(self):
        schema = {"type": "object", "properties": {"bounds": {"type": "array"}}}
        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema=schema,
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="test",
        )
        assert desc.show_options() == schema

    def test_validate_passes_through_dict(self):
        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="test",
        )
        params = {"bounds": [1.0, 2.0, 3.0, 4.0], "format": "xdmf"}
        result = desc.validate(params)
        assert result == params

    def test_validate_serializes_bounds_object(self):
        from dtcc_core.model import Bounds

        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="test",
        )
        bounds = Bounds(xmin=1.0, ymin=2.0, xmax=3.0, ymax=4.0)
        params = {"bounds": bounds}
        result = desc.validate(params)
        assert isinstance(result["bounds"], tuple)
        assert result["bounds"] == (1.0, 2.0, 3.0, 4.0)

    def test_call_injects_default_format(self):
        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf", "vtk"],
            source_service="test",
        )
        # Mock build to capture what it receives
        called_with = {}

        def mock_build(args, progress_callback=None, remote_info_callback=None):
            called_with.update(args)
            return b"data", "xdmf", "application/x-hdf5"

        desc.build = mock_build
        desc(bounds=[1.0, 2.0, 3.0, 4.0])
        assert called_with["format"] == "xdmf"

    def test_call_preserves_explicit_format(self):
        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf", "vtk"],
            source_service="test",
        )
        called_with = {}

        def mock_build(args, progress_callback=None, remote_info_callback=None):
            called_with.update(args)
            return b"data", "vtk", "application/x-vtk"

        desc.build = mock_build
        desc(bounds=[1.0, 2.0, 3.0, 4.0], format="vtk")
        assert called_with["format"] == "vtk"

    def test_args_model_accepts_any_kwargs(self):
        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="test",
        )
        # Should not raise -- _PassthroughArgs accepts anything
        instance = desc.ArgsModel(bounds=[1, 2, 3, 4], kappa=0.5, whatever="yes")
        assert instance.bounds == [1, 2, 3, 4]


def _mock_discovery_response():
    return {
        "service": "dtcc-sim",
        "version": "0.1.0",
        "protocol_version": "1",
        "datasets": {
            "mock_sim_dataset": {
                "name": "mock_sim_dataset",
                "description": "A mock simulation dataset",
                "args_schema": {"type": "object", "properties": {}},
                "result_kind": "mesh",
                "supported_formats": ["xdmf"],
                "timeout_hint": 600,
            }
        },
    }


class TestRegisterRemoteService:
    def test_success(self):
        from dtcc_core.datasets.remote import register_remote_service
        from dtcc_core.datasets import list as list_datasets, unregister

        mock_resp = MagicMock()
        mock_resp.json.return_value = _mock_discovery_response()
        mock_resp.raise_for_status = MagicMock()

        with patch("httpx.get", return_value=mock_resp):
            registered = register_remote_service("http://localhost:8001", timeout=5)

        assert "mock_sim_dataset" in registered
        assert "mock_sim_dataset" in list_datasets()

        ds = list_datasets()["mock_sim_dataset"]
        assert ds.source_service == "dtcc-sim"
        assert ds.data_category == "simulation"
        assert ds.timeout_hint == 600

        unregister("mock_sim_dataset")

    def test_unreachable_returns_empty(self):
        from dtcc_core.datasets.remote import register_remote_service
        import httpx

        with patch("httpx.get", side_effect=httpx.ConnectError("refused")):
            registered = register_remote_service("http://unreachable:9999", timeout=1)

        assert registered == []

    def test_500_error_returns_empty(self):
        from dtcc_core.datasets.remote import register_remote_service

        with patch("httpx.get", side_effect=Exception("500 Internal Server Error")):
            registered = register_remote_service("http://broken:8001", timeout=1)

        assert registered == []

    def test_cached_discoveries_roundtrip(self):
        from dtcc_core.datasets.remote import (
            register_remote_service,
            register_remote_descriptors_from_cache,
            get_cached_discoveries,
            _cached_service_discoveries,
        )
        from dtcc_core.datasets import list as list_datasets, unregister

        # Clear cache from prior test runs
        _cached_service_discoveries.clear()

        mock_resp = MagicMock()
        mock_resp.json.return_value = _mock_discovery_response()
        mock_resp.raise_for_status = MagicMock()

        with patch("httpx.get", return_value=mock_resp):
            register_remote_service("http://localhost:8001")

        cached = get_cached_discoveries()
        assert "http://localhost:8001" in cached

        # Unregister, then re-register from cache
        unregister("mock_sim_dataset")
        assert "mock_sim_dataset" not in list_datasets()

        register_remote_descriptors_from_cache(cached)
        assert "mock_sim_dataset" in list_datasets()

        # Cleanup
        unregister("mock_sim_dataset")
        _cached_service_discoveries.clear()


class TestRemoteValidationError:
    def test_422_raises_remote_validation_error(self):
        from dtcc_core.datasets.remote import RemoteValidationError

        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="test",
        )

        mock_resp = MagicMock()
        mock_resp.status_code = 422
        mock_resp.json.return_value = {"detail": [{"msg": "bounds required"}]}

        with patch("httpx.post", return_value=mock_resp):
            with pytest.raises(RemoteValidationError) as exc_info:
                desc.build({"bounds": [1, 2, 3, 4]})
            assert exc_info.value.status_code == 422
            assert "bounds required" in str(exc_info.value.detail)

    def test_remote_validation_error_is_value_error(self):
        """RemoteValidationError should be catchable as ValueError for compatibility."""
        from dtcc_core.datasets.remote import RemoteValidationError

        err = RemoteValidationError("test detail")
        assert isinstance(err, ValueError)


class TestPathContainment:
    def test_rejects_path_traversal(self):
        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="test",
        )

        # Mock submit to return a task_id
        mock_submit = MagicMock()
        mock_submit.status_code = 201
        mock_submit.json.return_value = {"task_id": "abc-123", "status": "pending"}
        mock_submit.raise_for_status = MagicMock()

        # Mock _stream_status to return a traversal path
        with patch("httpx.post", return_value=mock_submit):
            with patch.object(desc, "_stream_status", return_value="../../etc/passwd"):
                with pytest.raises(RuntimeError, match="unsafe result_file path"):
                    desc.build({"bounds": [1, 2, 3, 4]})

    def test_rejects_absolute_path(self):
        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="test",
        )

        mock_submit = MagicMock()
        mock_submit.status_code = 201
        mock_submit.json.return_value = {"task_id": "abc-123", "status": "pending"}
        mock_submit.raise_for_status = MagicMock()

        with patch("httpx.post", return_value=mock_submit):
            with patch.object(desc, "_stream_status", return_value="/etc/passwd"):
                with pytest.raises(RuntimeError, match="unsafe result_file path"):
                    desc.build({"bounds": [1, 2, 3, 4]})

    def test_accepts_safe_filename(self, tmp_path):
        import dtcc_core.datasets.remote as remote_mod

        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="test",
        )

        # Write a fake result file
        result_file = "abc-123.xdmf"
        fake_data = b"fake xdmf data"
        (tmp_path / result_file).write_bytes(fake_data)

        mock_submit = MagicMock()
        mock_submit.status_code = 201
        mock_submit.json.return_value = {"task_id": "abc-123", "status": "pending"}
        mock_submit.raise_for_status = MagicMock()

        with patch("httpx.post", return_value=mock_submit):
            with patch.object(desc, "_stream_status", return_value=result_file):
                with patch.object(remote_mod, "SHARED_RESULTS_DIR", str(tmp_path)):
                    data, ext, ct = desc.build({"bounds": [1, 2, 3, 4]})

        assert data == fake_data
        assert ext == "xdmf"
        assert not (tmp_path / result_file).exists()  # cleaned up


class TestRemoteInfoCallback:
    def test_remote_info_callback_called_with_task_id(self, tmp_path):
        import dtcc_core.datasets.remote as remote_mod

        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="test",
        )

        # Write a fake result file
        (tmp_path / "result.xdmf").write_bytes(b"fake data")

        mock_submit = MagicMock()
        mock_submit.status_code = 201
        mock_submit.json.return_value = {"task_id": "task-xyz", "status": "pending"}
        mock_submit.raise_for_status = MagicMock()

        callback = MagicMock()

        with patch("httpx.post", return_value=mock_submit):
            with patch.object(desc, "_stream_status", return_value="result.xdmf"):
                with patch.object(remote_mod, "SHARED_RESULTS_DIR", str(tmp_path)):
                    desc.build(
                        {"bounds": [1, 2, 3, 4]},
                        remote_info_callback=callback,
                    )

        callback.assert_called_once()
        info = callback.call_args[0][0]
        assert info["remote_task_id"] == "task-xyz"
        assert "cancel" in info["cancel_url"]


class TestSSEAndPolling:
    def test_poll_status_returns_result_file(self):
        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="test",
        )

        mock_resp = MagicMock()
        mock_resp.json.return_value = {
            "status": "completed",
            "result_file": "abc.xdmf",
            "size_bytes": 1024,
        }
        mock_resp.raise_for_status = MagicMock()

        with patch("httpx.get", return_value=mock_resp):
            result = desc._poll_status("http://localhost:8001/status/abc", None)

        assert result == "abc.xdmf"

    def test_poll_status_raises_on_failure(self):
        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="test",
        )

        mock_resp = MagicMock()
        mock_resp.json.return_value = {
            "status": "failed",
            "error": "Out of memory",
        }
        mock_resp.raise_for_status = MagicMock()

        with patch("httpx.get", return_value=mock_resp):
            with pytest.raises(RuntimeError, match="Out of memory"):
                desc._poll_status("http://localhost:8001/status/abc", None)

    def test_poll_status_forwards_progress(self):
        desc = RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="test",
        )

        call_count = 0
        progress_calls = []

        def mock_get(*args, **kwargs):
            nonlocal call_count
            call_count += 1
            resp = MagicMock()
            resp.raise_for_status = MagicMock()
            if call_count < 3:
                resp.json.return_value = {
                    "status": "running",
                    "progress": 0.5,
                    "message": "Building mesh",
                }
            else:
                resp.json.return_value = {
                    "status": "completed",
                    "result_file": "abc.xdmf",
                    "size_bytes": 1024,
                }
            return resp

        def progress_cb(state):
            progress_calls.append(state)

        with patch("httpx.get", side_effect=mock_get):
            with patch("time.sleep"):  # don't actually sleep
                result = desc._poll_status("http://localhost:8001/status/abc", progress_cb)

        assert result == "abc.xdmf"
        assert len(progress_calls) == 2
        assert progress_calls[0]["percent"] == 50.0  # 0.5 * 100


class TestStreamStatusFallback:
    def _make_desc(self):
        return RemoteDatasetDescriptor(
            name="test",
            description="",
            args_schema={},
            base_url="http://localhost:8001",
            result_kind="mesh",
            supported_formats=["xdmf"],
            source_service="test",
        )

    def test_sse_drop_falls_back_to_polling(self):
        """Unexpected SSE stream end should fall back to polling, not crash."""
        from dtcc_core.datasets.remote import SSEStreamError

        desc = self._make_desc()

        with patch.object(desc, "_stream_sse", side_effect=SSEStreamError("stream dropped")):
            with patch.object(desc, "_poll_status", return_value="abc.xdmf") as mock_poll:
                result = desc._stream_status("task-123", None)

        assert result == "abc.xdmf"
        mock_poll.assert_called_once()

    def test_transport_error_falls_back_to_polling(self):
        """httpx transport errors should fall back to polling."""
        import httpx

        desc = self._make_desc()

        with patch.object(desc, "_stream_sse", side_effect=httpx.ReadError("connection reset")):
            with patch.object(desc, "_poll_status", return_value="abc.xdmf") as mock_poll:
                result = desc._stream_status("task-123", None)

        assert result == "abc.xdmf"
        mock_poll.assert_called_once()

    def test_remote_protocol_error_falls_back_to_polling(self):
        """httpx.RemoteProtocolError (premature peer close) should fall back to polling."""
        import httpx

        desc = self._make_desc()

        with patch.object(desc, "_stream_sse", side_effect=httpx.RemoteProtocolError("peer closed")):
            with patch.object(desc, "_poll_status", return_value="abc.xdmf") as mock_poll:
                result = desc._stream_status("task-123", None)

        assert result == "abc.xdmf"
        mock_poll.assert_called_once()

    def test_terminal_failure_propagates_immediately(self):
        """RuntimeError from failed/cancelled jobs should NOT fall back to polling."""
        desc = self._make_desc()

        with patch.object(desc, "_stream_sse", side_effect=RuntimeError("Remote job failed: OOM")):
            with patch.object(desc, "_poll_status") as mock_poll:
                with pytest.raises(RuntimeError, match="OOM"):
                    desc._stream_status("task-123", None)

        mock_poll.assert_not_called()

    def test_terminal_cancel_propagates_immediately(self):
        """RuntimeError from cancellation should NOT fall back to polling."""
        desc = self._make_desc()

        with patch.object(desc, "_stream_sse", side_effect=RuntimeError("Remote job was cancelled")):
            with patch.object(desc, "_poll_status") as mock_poll:
                with pytest.raises(RuntimeError, match="cancelled"):
                    desc._stream_status("task-123", None)

        mock_poll.assert_not_called()
