"""Test Point geometry protobuf roundtrip."""

import pytest
import tempfile

# Import from dtcc_core package
from dtcc_core.model.geometry import Point
from dtcc_core.model.geometry.bounds import Bounds


def test_point_creation():
    """Test basic Point creation."""
    p = Point(x=1.5, y=2.5, z=3.5)
    assert p.x == 1.5
    assert p.y == 2.5
    assert p.z == 3.5


def test_point_bounds():
    """Test Point bounds calculation."""
    p = Point(x=10.0, y=20.0, z=30.0)
    b = p.bounds
    
    assert b.xmin == 10.0
    assert b.xmax == 10.0
    assert b.ymin == 20.0
    assert b.ymax == 20.0
    assert b.zmin == 30.0
    assert b.zmax == 30.0


def test_point_offset():
    """Test Point offset method."""
    p = Point(x=1.0, y=2.0, z=3.0)
    p.offset(dx=1.0, dy=2.0, dz=3.0)
    
    assert p.x == 2.0
    assert p.y == 4.0
    assert p.z == 6.0


def test_point_protobuf_roundtrip():
    """Test Point protobuf serialization and deserialization."""
    # Create original point
    p1 = Point(x=12.34, y=56.78, z=90.12)
    
    # Serialize to protobuf
    pb = p1.to_proto()
    
    # Deserialize back
    p2 = Point()
    p2.from_proto(pb)
    
    # Verify coordinates match
    assert p2.x == pytest.approx(12.34)
    assert p2.y == pytest.approx(56.78)
    assert p2.z == pytest.approx(90.12)


def test_point_protobuf_bytes_roundtrip():
    """Test Point protobuf roundtrip through bytes."""
    # Create original point
    p1 = Point(x=111.111, y=222.222, z=333.333)
    
    # Serialize to protobuf and then to bytes
    pb = p1.to_proto()
    pb_bytes = pb.SerializeToString()
    
    # Deserialize from bytes
    p2 = Point()
    p2.from_proto(pb_bytes)
    
    # Verify coordinates match
    assert p2.x == pytest.approx(111.111)
    assert p2.y == pytest.approx(222.222)
    assert p2.z == pytest.approx(333.333)


def test_point_default_values():
    """Test Point with default z value."""
    p = Point(x=1.0, y=2.0)  # z defaults to 0
    
    assert p.x == 1.0
    assert p.y == 2.0
    assert p.z == 0.0
    
    # Test roundtrip preserves default
    pb = p.to_proto()
    p2 = Point()
    p2.from_proto(pb)
    
    assert p2.z == 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
