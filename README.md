This is library written in C++ which contains some useful structures and functions that are based on linear algebra and are commonly used in computer graphics.	
Library provides an easy way to do vector dot products, cross products, vector-matrix multiplications and many more.
Structures: 
	vector2, vector3, vector4, matrix2, matrix3, matrix4.
Functions:
	//Function creates rotation matrix to rotate 3D vector along the Z axis at an angle "angle".
	matrix4 RotationMatrixZ(float angle);

	//Function creates rotation matrix to rotate 3D vector along the Y axis at an angle "angle".
	matrix4 RotationMatrixY(float angle);

	//Function creates rotation matrix to rotate 3D vector along the X axis at an angle "angle".
	matrix4 RotationMatrixX(float angle);

	//Function creates rotation matrix to rotate 3D vector along the "axis" axis at an angle "angle".
	matrix4 RotationMatrixAxis(const vector3& origin, const vector3& axis, float angle);

	//Function creates rotation matrix to redirect the 3D vector along the "newDir".
	matrix4 RotationMatrixRedirect(vector3 oldDir, vector3 newDir);

	//Function creates matrix to scale 3D vector along x, y, z axes.
	matrix4 ScaleMatrix(float x, float y, float z);

	//Function creates matrix to scale 3D vector along "Axis" axis.
	matrix4 ScaleMatrixAxis(vector3 axis, float k);

	//Function creates matrix to offset 3D vector along x, y, z axes.
	matrix4 OffsetMatrix(float offsetX, float offsetY, float offsetZ);

	//Function creates matrix to translate 3D vector to camera's view space.
	matrix4 LookToMatrix(const vector3 & pos, const vector3 & dir, const vector3 & up);

	//Function creates matrix to translate 3D vector to camera's view space.
	matrix4 LookAtMatrix(const vector3 & pos, const vector3 & at, const vector3 & up);

	//Function creates perspective projection matrix. It projects vector on the "nearplane".
	matrix4 ProjectionMatrixPerspective(float fov, float nearPlane, float farPlane, float aspectRatio);

	//Function creates orthogonal projection matrix.
	matrix4 ProjectionMatrixOrthogonal(float width, float height, float farPlane);

	//Function creates matrix, multiplying by which reflects 3D vector from the plane defined by its "normal" vector.
	matrix4 ReflectMatrix(vector3 normal);

	//Function generates tangent and bitangent vectors for each vertex in the array based on its texture coordinates. Vertex array mustn't be indexed!
	bool GenerateVertexTangents(const vector3* vertices, const vector2* texCoords, vector3* tangents, vector3* bitangents, int count);

	//Function generates tangent and bitangent vectors for each vertex in the array based on its texture coordinates and position. Vertex array must be indexed!
	bool GenerateVertexTangentsIndexed(const vector3* vertices, const vector2* texCoords, vector3* tangents, vector3* bitangents, unsigned count, const unsigned* indices, unsigned indicesCount);

	//Function generates normal for each 3 vertices in the array based on their position.
	bool GenerateVertexNormalsCCW(const vector3* vertices, vector3* normals, unsigned count);

	//Function generates normal for each 3 vertices in the array based on their position.
	bool GenerateVertexNormalsCW(const vector3* vertices, vector3* normals, unsigned count);

	//Function generates normal for each vertex in the array based on its position.
	bool GenerateVertexNormalsCCWIndexed(const vector3* vertices, vector3* normals, unsigned count, const unsigned* indices, unsigned indicesCount);

	//Function generates normal for each vertex in the array based on its position.
	bool GenerateVertexNormalsCWIndexed(const vector3* vertices, vector3* normals, unsigned count, const unsigned* indices, unsigned indicesCount);

	//Function generates texture coordinates for mesh. Texture is projected spherically.
	void GenerateSphericalTexCoords(const vector3* vertices, vector2* coords, unsigned count, float uOffset, float vOffset, bool uInverse, float uTile, float vTile);
	