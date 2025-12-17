from dtcc_core.model import Mesh

import numpy as np
import numpy as np


def create_sphere(center, radius, resolution=20, scale=(1.0, 1.0, 1.0)):
    """
    Create a triangulated sphere mesh with optional scaling along each axis.

    Parameters:
    -----------
    center : list or numpy.ndarray
        The (x, y, z) coordinates of the sphere center.
    radius : float
        The base radius of the sphere before scaling.
    resolution : int, optional
        The resolution of the sphere mesh. Higher values create more triangles.
    scale : tuple of 3 floats, optional
        Scaling factors for the x, y, and z axes respectively.
        Use values < 1 to squash, values > 1 to stretch.
        Default is (1.0, 1.0, 1.0) for a perfect sphere.

    Returns:
    --------
    vertices : numpy.ndarray
        An Nx3 array of vertex coordinates.
    faces : numpy.ndarray
        An Mx3 array of face indices, where each row contains the indices of
        the vertices that make up a triangular face.
    """
    # Convert center to numpy array if it isn't already
    center = np.array(center)
    scale = np.array(scale)

    # Create a grid of points on a unit sphere
    phi = np.linspace(0, 2 * np.pi, resolution)
    theta = np.linspace(0, np.pi, resolution)

    # Create the mesh grid
    phi, theta = np.meshgrid(phi, theta)

    # Convert to Cartesian coordinates on a unit sphere
    x = np.cos(phi) * np.sin(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(theta)

    # Reshape to get a list of 3D points
    x = x.flatten()
    y = y.flatten()
    z = z.flatten()

    # Stack coordinates
    vertices = np.column_stack((x, y, z))

    # Apply scaling factors to each axis
    vertices = vertices * scale

    # Apply radius and center offset
    vertices = vertices * radius + center

    # Create faces
    num_points = resolution * resolution
    faces = []

    for i in range(resolution):
        for j in range(resolution):
            p1 = i * resolution + j
            p2 = (i + 1) % resolution * resolution + j
            p3 = (i + 1) % resolution * resolution + (j + 1) % resolution
            p4 = i * resolution + (j + 1) % resolution

            # Add two triangular faces for each grid cell
            if j != resolution - 1 and i != resolution - 1:
                faces.append([p1, p2, p3])
                faces.append([p1, p3, p4])

    mesh = Mesh(vertices=vertices, faces=np.array(faces))

    return mesh


import numpy as np


def create_cylinder(center, radius, height, resolution=20, axis=2, cap=True):
    """
    Create a triangulated cylinder mesh.

    Parameters:
    -----------
    center : list or numpy.ndarray
        The (x, y, z) coordinates of the cylinder center.
    radius : float
        The radius of the cylinder.
    height : float
        The height of the cylinder.
    resolution : int, optional
        The resolution of the cylinder mesh (number of points around circumference).
        Higher values create more triangles.
    axis : int, optional
        The axis along which to align the cylinder's height:
        0 for x-axis, 1 for y-axis, 2 for z-axis (default).
    cap : bool, optional
        Whether to include the end caps of the cylinder.

    Returns:
    --------
    vertices : numpy.ndarray
        An Nx3 array of vertex coordinates.
    faces : numpy.ndarray
        An Mx3 array of face indices, where each row contains the indices of
        the vertices that make up a triangular face.
    """
    # Convert center to numpy array if it isn't already
    center = np.array(center)

    # Generate circle points on the unit circle
    theta = np.linspace(0, 2 * np.pi, resolution, endpoint=False)
    x = np.cos(theta)
    y = np.sin(theta)

    # Initialize vertices and faces lists
    vertices = []
    faces = []

    # Create the cylinder sides
    for i in range(resolution):
        # For each angle, create points at the top and bottom of cylinder
        for h in [-0.5, 0.5]:  # Normalized height positions
            # Create the correct coordinate based on the axis parameter
            if axis == 0:  # X-axis aligned
                point = np.array([h * height, x[i] * radius, y[i] * radius])
            elif axis == 1:  # Y-axis aligned
                point = np.array([x[i] * radius, h * height, y[i] * radius])
            else:  # Z-axis aligned (default)
                point = np.array([x[i] * radius, y[i] * radius, h * height])

            vertices.append(center + point)

    # Create side faces (rectangles split into triangles)
    for i in range(resolution):
        # Define the four corners of each rectangle
        i0 = i * 2  # Bottom of current segment
        i1 = (i * 2) + 1  # Top of current segment
        i2 = ((i + 1) % resolution) * 2  # Bottom of next segment
        i3 = ((i + 1) % resolution) * 2 + 1  # Top of next segment

        # Add two triangular faces
        faces.append([i0, i2, i3])
        faces.append([i0, i3, i1])

    # Add caps if requested
    if cap:
        # Add center points for top and bottom caps
        bottom_center_idx = len(vertices)
        if axis == 0:
            bottom_center = center + np.array([-0.5 * height, 0, 0])
            top_center = center + np.array([0.5 * height, 0, 0])
        elif axis == 1:
            bottom_center = center + np.array([0, -0.5 * height, 0])
            top_center = center + np.array([0, 0.5 * height, 0])
        else:  # axis == 2
            bottom_center = center + np.array([0, 0, -0.5 * height])
            top_center = center + np.array([0, 0, 0.5 * height])

        vertices.append(bottom_center)
        top_center_idx = len(vertices)
        vertices.append(top_center)

        # Add triangular faces for bottom cap
        for i in range(resolution):
            i0 = i * 2  # Current bottom point
            i1 = ((i + 1) % resolution) * 2  # Next bottom point
            faces.append([bottom_center_idx, i1, i0])

        # Add triangular faces for top cap
        for i in range(resolution):
            i0 = i * 2 + 1  # Current top point
            i1 = ((i + 1) % resolution) * 2 + 1  # Next top point
            faces.append([top_center_idx, i0, i1])

    mesh = Mesh(vertices=np.array(vertices), faces=np.array(faces))

    return mesh


import numpy as np


def create_cone_mesh(center, radius, height, resolution=20, axis=2, cap=True):
    """
    Create a triangulated cone mesh.

    Parameters:
    -----------
    center : list or numpy.ndarray
        The (x, y, z) coordinates of the cone's base center.
    radius : float
        The radius of the cone base.
    height : float
        The height of the cone from base to tip.
    resolution : int, optional
        The resolution of the cone mesh (number of points around base circumference).
        Higher values create more triangles.
    axis : int, optional
        The axis along which to align the cone's height:
        0 for x-axis, 1 for y-axis, 2 for z-axis (default).
    cap : bool, optional
        Whether to include the base cap of the cone.

    Returns:
    --------
    vertices : numpy.ndarray
        An Nx3 array of vertex coordinates.
    faces : numpy.ndarray
        An Mx3 array of face indices, where each row contains the indices of
        the vertices that make up a triangular face.
    """
    # Convert center to numpy array if it isn't already
    center = np.array(center)

    # Generate circle points for the base
    theta = np.linspace(0, 2 * np.pi, resolution, endpoint=False)
    x = np.cos(theta) * radius
    y = np.sin(theta) * radius

    # Initialize vertices and faces lists
    vertices = []
    faces = []

    # Create base vertices
    for i in range(resolution):
        # Create base points based on the axis parameter
        if axis == 0:  # X-axis aligned
            point = np.array([0, x[i], y[i]])
        elif axis == 1:  # Y-axis aligned
            point = np.array([x[i], 0, y[i]])
        else:  # Z-axis aligned (default)
            point = np.array([x[i], y[i], 0])

        vertices.append(center + point)

    # Determine the tip position
    if axis == 0:
        tip = center + np.array([height, 0, 0])
    elif axis == 1:
        tip = center + np.array([0, height, 0])
    else:  # axis == 2
        tip = center + np.array([0, 0, height])

    # Add tip vertex
    tip_idx = len(vertices)
    vertices.append(tip)

    # Create triangular faces for the cone surface
    for i in range(resolution):
        i0 = i  # Current base point
        i1 = (i + 1) % resolution  # Next base point

        # Add triangular face from base edge to tip
        faces.append([i0, i1, tip_idx])

    # Add base cap if requested
    if cap:
        # Add center point for base cap
        base_center_idx = len(vertices)
        vertices.append(center)

        # Add triangular faces for base cap
        for i in range(resolution):
            i0 = i  # Current base point
            i1 = (i + 1) % resolution  # Next base point
            faces.append([base_center_idx, i0, i1])

    mesh = Mesh(vertices=np.array(vertices), faces=np.array(faces))
    return mesh
