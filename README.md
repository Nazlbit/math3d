This is a library written in C++ that contains some useful structures and functions that are based on linear algebra and are commonly used in computer graphics.	

The library provides an easy way to do vector dot products, cross products, vector-matrix multiplications and more.

Structures: vector2, vector3, vector4, matrix2, matrix3, matrix4.

Functions:
	
	//The function creates a rotation matrix to rotate a 3D vector along the Z axis at an angle "angle".
	matrix4 RotationMatrixZ(float angle);
	
	//The function creates a rotation matrix to rotate a 3D vector along the Y axis at an angle "angle".
	matrix4 RotationMatrixY(float angle);

	//The function creates a rotation matrix to rotate a 3D vector along the X axis at an angle "angle".
	matrix4 RotationMatrixX(float angle);

	//The function creates a rotation matrix to rotate a 3D vector along the "axis" axis at an angle "angle".
	matrix4 RotationMatrixAxis(const vector3& origin, const vector3& axis, float angle);

	//The function creates a rotation matrix to redirect the 3D vector along the "newDir".
	matrix4 RotationMatrixRedirect(vector3 oldDir, vector3 newDir);

	//The function creates a matrix to scale a 3D vector along the x, y and z axes.
	matrix4 ScaleMatrix(float x, float y, float z);

	//The function creates a matrix to scale a 3D vector along the "Axis" axis.
	matrix4 ScaleMatrixAxis(vector3 axis, float k);

	//The function creates a matrix to offset a 3D vector along the x, y and z axes.
	matrix4 OffsetMatrix(float offsetX, float offsetY, float offsetZ);

	//The function creates a matrix to translate a 3D vector to the view space of a virtual camera.
	matrix4 LookToMatrix(const vector3 & pos, const vector3 & dir, const vector3 & up);

	//The function creates a matrix to translate a 3D vector to the view space of a virtual camera.
	matrix4 LookAtMatrix(const vector3 & pos, const vector3 & at, const vector3 & up);

	//The function creates a perspective projection matrix.
	matrix4 ProjectionMatrixPerspective(float fov, float nearPlane, float farPlane, float aspectRatio);

	//The function creates an orthogonal projection matrix.
	matrix4 ProjectionMatrixOrthogonal(float width, float height, float farPlane);

	//The function creates a matrix, multiplying by which reflects the 3D vector from the plane defined by its "normal" vector.
	matrix4 ReflectMatrix(vector3 normal);

	//The function generates the tangent and bitangent vectors for each vertex in the array based on its texture coordinates. The vertex array mustn't be indexed!
	bool GenerateVertexTangents(const vector3* vertices, const vector2* texCoords, vector3* tangents, vector3* bitangents, int count);

	//The function generates the tangent and bitangent vectors for each vertex in the array based on its texture coordinates and position. The vertex array must be indexed!
	bool GenerateVertexTangentsIndexed(const vector3* vertices, const vector2* texCoords, vector3* tangents, vector3* bitangents, unsigned count, const unsigned* indices, unsigned indicesCount);

	//The function generates normals for the triangles defined by the vertices in the array based on their position.
	bool GenerateVertexNormalsCCW(const vector3* vertices, vector3* normals, unsigned count);

	//The function generates normals for the triangles defined by the vertices in the array based on their position.
	bool GenerateVertexNormalsCW(const vector3* vertices, vector3* normals, unsigned count);

	//The function generates normals for the triangles defined by the vertices in the array based on their position.
	bool GenerateVertexNormalsCCWIndexed(const vector3* vertices, vector3* normals, unsigned count, const unsigned* indices, unsigned indicesCount);

	//The function generates normals for the triangles defined by the vertices in the array based on their position.
	bool GenerateVertexNormalsCWIndexed(const vector3* vertices, vector3* normals, unsigned count, const unsigned* indices, unsigned indicesCount);

	//The function generates texture coordinates for the mesh. The texture is projected spherically.
	void GenerateSphericalTexCoords(const vector3* vertices, vector2* coords, unsigned count, float uOffset, float vOffset, bool uInverse, float uTile, float vTile);
	
