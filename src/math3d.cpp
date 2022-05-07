#include "math3d.h"
#include <cmath>

namespace math3d
{

	vector2::vector2() {}

	vector2::vector2(float _x, float _y)
	{
		x = _x;
		y = _y;
	}

	vector2::vector2(const vector3 &a)
	{
		x = a.x;
		y = a.y;
	}

	vector2::vector2(const vector4 &a)
	{
		x = a.x;
		y = a.y;
	}

	vector2 vector2::operator+(const vector2 &a) const
	{
		return vector2{x + a.x, y + a.y};
	}

	bool vector2::operator==(const vector2 &a) const
	{
		if (x == a.x && y == a.y)
			return true;
		return false;
	}

	vector2 &vector2::operator+=(const vector2 &a)
	{
		x += a.x;
		y += a.y;

		return *this;
	}

	vector2 &vector2::operator-=(const vector2 &a)
	{
		x -= a.x;
		y -= a.y;

		return *this;
	}

	vector2 vector2::operator-(const vector2 &a) const
	{
		return vector2{x - a.x, y - a.y};
	}

	vector2 vector2::operator-() const
	{
		return vector2{-x, -y};
	}

	vector2 vector2::operator*(float a) const
	{
		return vector2{x * a, y * a};
	}

	vector2::operator float *()
	{
		return all;
	}

	float vector2::operator*(const vector2 &a) const
	{
		return x * a.x + y * a.y;
	}

	vector2 vector2::operator/(float a) const
	{
		return vector2{x / a, y / a};
	}

	float vector2::Length() const
	{
		return sqrtf(x * x + y * y);
	}

	vector2 vector2::Normalize() const
	{
		return *this / Length();
	}

	vector3::vector3() {}

	vector3::vector3(float _x, float _y, float _z)
	{
		x = _x;
		y = _y;
		z = _z;
	}

	vector3::vector3(const vector2 &_xy, float _z)
	{
		x = _xy.x;
		y = _xy.y;
		z = _z;
	}

	vector3::vector3(const vector4 &_xyz)
	{
		x = _xyz.x;
		y = _xyz.y;
		z = _xyz.z;
	}

	vector3 vector3::operator+(const vector3 &a) const
	{
		return vector3{x + a.x, y + a.y, z + a.z};
	}

	vector3 vector3::operator-(const vector3 &a) const
	{
		return vector3{x - a.x, y - a.y, z - a.z};
	}

	vector3 vector3::operator-() const
	{
		return vector3{-x, -y, -z};
	}

	vector3 vector3::operator*(float a) const
	{
		return vector3{x * a, y * a, z * a};
	}

	vector3::operator float *()
	{
		return all;
	}

	float vector3::operator*(const vector3 &a) const
	{
		return x * a.x + y * a.y + z * a.z;
	}

	vector3 vector3::operator^(const vector3 &a) const
	{
		return vector3{y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x};
	}

	vector3 vector3::operator/(float a) const
	{
		return vector3{x / a, y / a, z / a};
	}

	float vector3::Length() const
	{
		return sqrtf(x * x + y * y + z * z);
	}

	vector3 vector3::Normalize() const
	{
		return *this / Length();
	}

	bool vector3::operator==(const vector3 &a) const
	{
		if (x == a.x && y == a.y && z == a.z)
			return true;
		return false;
	}

	vector3 &vector3::operator+=(const vector3 &a)
	{
		x += a.x;
		y += a.y;
		z += a.z;

		return *this;
	}

	vector3 &vector3::operator-=(const vector3 &a)
	{
		x -= a.x;
		y -= a.y;
		z -= a.z;

		return *this;
	}

	vector4::vector4() {}

	vector4::vector4(float _x, float _y, float _z, float _w)
	{
		x = _x;
		y = _y;
		z = _z;
		w = _w;
	}

	vector4::vector4(vector3 _xyz, float _w)
	{
		x = _xyz.x;
		y = _xyz.y;
		z = _xyz.z;
		w = _w;
	}

	vector4::vector4(vector2 _xy, float _z, float _w)
	{
		x = _xy.x;
		y = _xy.y;
		z = _z;
		w = _w;
	}

	vector4 vector4::operator+(const vector4 &a) const
	{
		return vector4{x + a.x, y + a.y, z + a.z, w + a.w};
	}

	vector4 vector4::operator-(const vector4 &a) const
	{
		return vector4{x - a.x, y - a.y, z - a.z, w - a.w};
	}

	vector4 vector4::operator-() const
	{
		return vector4{-x, -y, -z, -w};
	}

	vector4 vector4::operator*(float a) const
	{
		return vector4{x * a, y * a, z * a, w * a};
	}

	vector4::operator float *()
	{
		return all;
	}

	float vector4::operator*(const vector4 &a) const
	{
		return x * a.x + y * a.y + z * a.z + w * a.w;
	}

	vector4 vector4::operator/(float a) const
	{
		return vector4{x / a, y / a, z / a, w / a};
	}

	float vector4::Length() const
	{
		return sqrtf(x * x + y * y + z * z + w * w);
	}

	vector4 vector4::Normalize() const
	{
		return *this / Length();
	}

	bool vector4::operator==(const vector4 &a) const
	{
		if (x == a.x && y == a.y && z == a.z && w == a.w)
			return true;
		return false;
	}

	vector4 &vector4::operator+=(const vector4 &a)
	{
		x += a.x;
		y += a.y;
		z += a.z;
		w += a.w;

		return *this;
	}

	vector4 &vector4::operator-=(const vector4 &a)
	{
		x -= a.x;
		y -= a.y;
		z -= a.z;
		w -= a.w;

		return *this;
	}

	matrix2::matrix2() {}

	matrix2::matrix2(float _m11, float _m21, float _m12, float _m22)
	{
		m11 = _m11;
		m21 = _m21;
		m12 = _m12;
		m22 = _m22;
	}

	matrix2::matrix2(const vector2 &column1, const vector2 &column2)
	{
		columns[0] = column1;
		columns[1] = column2;
	}

	matrix2::matrix2(const matrix3 &a, unsigned offsetX, unsigned offsetY)
	{
		all2D[0][0] = a.all2D[offsetX][offsetY];
		all2D[0][1] = a.all2D[offsetX][offsetY + 1];
		all2D[1][0] = a.all2D[offsetX + 1][offsetY];
		all2D[1][1] = a.all2D[offsetX + 1][offsetY + 1];
	}

	matrix2::matrix2(const matrix4 &a, unsigned offsetX, unsigned offsetY)
	{
		all2D[0][0] = a.all2D[offsetX][offsetY];
		all2D[0][1] = a.all2D[offsetX][offsetY + 1];
		all2D[1][0] = a.all2D[offsetX + 1][offsetY];
		all2D[1][1] = a.all2D[offsetX + 1][offsetY + 1];
	}

	float matrix2::Determinant() const
	{
		return m11 * m22 - m21 * m12;
	}

	matrix2 matrix2::Transpose() const
	{
		return matrix2{m11, m12, m21, m22};
	}

	matrix2 matrix2::Inverse() const
	{
		float det = Determinant();
		if (!det)
			return *this;
		return matrix2{m22, -m21, -m12, m11} * (1.0f / det);
	}

	matrix2 matrix2::operator*(float a) const
	{
		return matrix2{m11 * a, m21 * a, m12 * a, m22 * a};
	}

	matrix2 matrix2::operator+(const matrix2 &a) const
	{
		return matrix2{m11 + a.m11, m21 + a.m21, m12 + a.m12, m22 + a.m22};
	}

	matrix2 matrix2::operator-(const matrix2 &a) const
	{
		return matrix2{m11 - a.m11, m21 - a.m21, m12 - a.m12, m22 - a.m22};
	}

	matrix2 matrix2::operator*(const matrix2 &a) const
	{
		matrix2 res;
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				res.all2D[j][i] = 0;
				for (int k = 0; k < 2; k++)
				{
					res.all2D[j][i] += all2D[k][i] * a.all2D[j][k];
				}
			}
		}
		return res;
	}

	vector2 &matrix2::operator[](int index)
	{
		return columns[index];
	}

	matrix3::matrix3() {}

	matrix3::matrix3(float _m11, float _m21, float _m31, float _m12, float _m22, float _m32, float _m13, float _m23, float _m33)
	{
		m11 = _m11;
		m21 = _m21;
		m31 = _m31;
		m12 = _m12;
		m22 = _m22;
		m32 = _m32;
		m13 = _m13;
		m23 = _m23;
		m33 = _m33;
	}

	matrix3::matrix3(const vector3 &column1, const vector3 &column2, const vector3 &column3)
	{
		columns[0] = column1;
		columns[1] = column2;
		columns[2] = column3;
	}

	matrix3::matrix3(const matrix2 &a, unsigned offsetX, unsigned offsetY)
	{
		all2D[offsetX][offsetY] = a.all2D[0][0];
		all2D[offsetX][offsetY + 1] = a.all2D[0][1];
		all2D[offsetX + 1][offsetY] = a.all2D[1][0];
		all2D[offsetX + 1][offsetY + 1] = a.all2D[1][1];
		all2D[2][2] = 1;
	}

	matrix3::matrix3(const matrix4 &a, unsigned offsetX, unsigned offsetY)
	{
		for (unsigned i = 0; i < 3; i++)
		{
			for (unsigned j = 0; j < 3; j++)
			{
				all2D[i][j] = a.all2D[i + offsetX][j + offsetY];
			}
		}
	}

	float matrix3::Determinant() const
	{
		return m11 * matrix2{m22, m32, m23, m33}.Determinant() - m12 * matrix2{m21, m31, m23, m33}.Determinant() + m13 * matrix2{m21, m31, m22, m32}.Determinant();
	}

	matrix3 matrix3::Transpose() const
	{
		return matrix3{m11, m12, m13, m21, m22, m23, m31, m32, m33};
	}

	matrix3 matrix3::operator*(float a) const
	{
		return matrix3{m11 * a, m21 * a, m31 * a, m12 * a, m22 * a, m32 * a, m13 * a, m23 * a, m33 * a};
	}

	matrix3 matrix3::operator+(const matrix3 &b) const
	{
		return matrix3{m11 + b.m11, m21 + b.m21, m31 + b.m31, m12 + b.m12, m22 + b.m22, m32 + b.m32, m13 + b.m13, m23 + b.m23, m33 + b.m33};
	}

	matrix3 matrix3::operator-(const matrix3 &b) const
	{
		return matrix3{m11 - b.m11, m21 - b.m21, m31 - b.m31, m12 - b.m12, m22 - b.m22, m32 - b.m32, m13 - b.m13, m23 - b.m23, m33 - b.m33};
	}

	matrix3 matrix3::operator*(const matrix3 &b) const
	{
		matrix3 res;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				res.all2D[j][i] = 0;
				for (int k = 0; k < 3; k++)
				{
					res.all2D[j][i] += all2D[k][i] * b.all2D[j][k];
				}
			}
		}
		return res;
	}

	matrix3 matrix3::Inverse() const
	{
		float det = Determinant();
		if (!det)
			return *this;
		return matrix3{matrix2{m22, m32, m23, m33}.Determinant(), -matrix2{m21, m31, m23, m33}.Determinant(), matrix2{m21, m31, m22, m32}.Determinant(),
					   -matrix2{m12, m32, m13, m33}.Determinant(), matrix2{m11, m31, m13, m33}.Determinant(), -matrix2{m11, m31, m12, m32}.Determinant(),
					   matrix2{m12, m22, m13, m23}.Determinant(), -matrix2{m11, m21, m13, m23}.Determinant(), matrix2{m11, m21, m12, m22}.Determinant()} *
			   (1.0f / det);
	}

	matrix3::operator float *()
	{
		return all;
	}

	vector3 &matrix3::operator[](int index)
	{
		return columns[index];
	}

	matrix4::matrix4() {}

	matrix4::matrix4(float _m11, float _m21, float _m31, float _m41, float _m12, float _m22, float _m32, float _m42, float _m13, float _m23, float _m33, float _m43, float _m14, float _m24, float _m34, float _m44)
	{
		m11 = _m11;
		m21 = _m21;
		m31 = _m31;
		m41 = _m41;
		m12 = _m12;
		m22 = _m22;
		m32 = _m32;
		m42 = _m42;
		m13 = _m13;
		m23 = _m23;
		m33 = _m33;
		m43 = _m43;
		m14 = _m14;
		m24 = _m24;
		m34 = _m34;
		m44 = _m44;
	}

	matrix4::matrix4(const matrix2 &a, unsigned offsetX, unsigned offsetY)
	{
		all2D[offsetX][offsetY] = a.all2D[0][0];
		all2D[offsetX][offsetY + 1] = a.all2D[0][1];
		all2D[offsetX + 1][offsetY] = a.all2D[1][0];
		all2D[offsetX + 1][offsetY + 1] = a.all2D[1][1];
		all2D[2][2] = 1;
		all2D[3][3] = 1;
	}

	matrix4::matrix4(const matrix3 &a, unsigned offsetX, unsigned offsetY)
	{
		for (unsigned i = 0; i < 3; i++)
		{
			for (unsigned j = 0; j < 3; j++)
			{
				all2D[i + offsetX][j + offsetY] = a.all2D[i][j];
			}
		}
		all2D[3][3] = 1;
	}

	matrix4::matrix4(const vector4 &column1, const vector4 &column2, const vector4 &column3, const vector4 &column4)
	{
		columns[0] = column1;
		columns[1] = column2;
		columns[2] = column3;
		columns[3] = column4;
	}

	float matrix4::Determinant() const
	{
		return m11 * matrix3{m22, m32, m42, m23, m33, m43, m24, m34, m44}.Determinant() -
			   m12 * matrix3{m21, m31, m41, m23, m33, m43, m24, m34, m44}.Determinant() +
			   m13 * matrix3{m21, m31, m41, m22, m32, m42, m24, m34, m44}.Determinant() -
			   m14 * matrix3{m21, m31, m41, m22, m32, m42, m23, m33, m43}.Determinant();
	}

	matrix4 matrix4::Transpose() const
	{
		return matrix4{m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44};
	}

	matrix4 matrix4::Inverse() const
	{
		float det = Determinant();
		if (!det)
			return *this;
		return matrix4{matrix3{m22, m32, m42, m23, m33, m43, m24, m34, m44}.Determinant(),
					   matrix3{m21, m31, m41, m23, m33, m43, m24, m34, m44}.Determinant(),
					   matrix3{m21, m31, m41, m22, m32, m42, m24, m34, m44}.Determinant(),
					   matrix3{m21, m31, m41, m22, m32, m42, m23, m33, m43}.Determinant(),

					   matrix3{m12, m32, m44, m13, m33, m43, m14, m34, m44}.Determinant(),
					   matrix3{m11, m31, m41, m13, m33, m43, m14, m34, m44}.Determinant(),
					   matrix3{m11, m31, m41, m12, m32, m42, m14, m34, m44}.Determinant(),
					   matrix3{m11, m31, m41, m12, m32, m42, m13, m33, m43}.Determinant(),

					   matrix3{m12, m22, m42, m13, m23, m43, m14, m24, m44}.Determinant(),
					   matrix3{m11, m21, m41, m13, m23, m43, m14, m24, m44}.Determinant(),
					   matrix3{m11, m21, m41, m12, m22, m42, m14, m24, m44}.Determinant(),
					   matrix3{m11, m21, m41, m12, m22, m42, m13, m23, m43}.Determinant(),

					   matrix3{m12, m22, m32, m13, m23, m33, m14, m24, m34}.Determinant(),
					   matrix3{m11, m21, m31, m13, m23, m33, m14, m24, m34}.Determinant(),
					   matrix3{m11, m21, m31, m12, m22, m32, m14, m24, m34}.Determinant(),
					   matrix3{m11, m21, m31, m12, m22, m32, m13, m23, m33}.Determinant()} *
			   (1.0f / det);
	}

	matrix4 matrix4::operator*(float a) const
	{
		return matrix4{a * m11, a * m21, a * m31, a * m41, a * m12, a * m22, a * m32, a * m42, a * m13, a * m23, a * m33, a * m43, a * m14, a * m24, a * m34, a * m44};
	}

	matrix4 matrix4::operator+(const matrix4 &b) const
	{
		matrix4 res = {m11 + b.m11, m21 + b.m21, m31 + b.m31, m41 + b.m41,
					   m12 + b.m12, m22 + b.m22, m32 + b.m32, m42 + b.m42,
					   m13 + b.m13, m23 + b.m23, m33 + b.m33, m43 + b.m43,
					   m14 + b.m14, m24 + b.m24, m34 + b.m34, m44 + b.m44};
		return res;
	}

	matrix4 matrix4::operator-(const matrix4 &b) const
	{
		matrix4 res = {m11 - b.m11, m21 - b.m21, m31 - b.m31, m41 - b.m41,
					   m12 - b.m12, m22 - b.m22, m32 - b.m32, m42 - b.m42,
					   m13 - b.m13, m23 - b.m23, m33 - b.m33, m43 - b.m43,
					   m14 - b.m14, m24 - b.m24, m34 - b.m34, m44 - b.m44};
		return res;
	}

	matrix4 matrix4::operator*(const matrix4 &b) const
	{
		matrix4 res;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				res.all2D[j][i] = 0;
				for (int k = 0; k < 4; k++)
				{
					res.all2D[j][i] += all2D[k][i] * b.all2D[j][k];
				}
			}
		}
		return res;
	}

	matrix4::operator float *()
	{
		return all;
	}

	vector4 &matrix4::operator[](int index)
	{
		return columns[index];
	}

	matrix4 operator*(float a, const matrix4 &b)
	{
		return b * a;
	}

	matrix3 operator*(float a, const matrix3 &b)
	{
		return b * a;
	}

	matrix2 operator*(float a, const matrix2 &b)
	{
		return b * a;
	}

	vector3 operator*(float a, const vector3 &b)
	{
		return b * a;
	}

	vector2 operator*(const vector2 &a, const matrix2 &b)
	{
		return vector2(a.x * b.m11 + a.y * b.m21,
					   a.x * b.m12 + a.y * b.m22);
	}

	vector2 operator*(const matrix2 &a, const vector2 &b)
	{
		return vector2(b.x * a.m11 + b.y * a.m12,
					   b.x * a.m21 + b.y * a.m22);
	}

	vector3 operator*(const vector3 &a, const matrix3 &b)
	{
		return vector3{a.x * b.m11 + a.y * b.m21 + a.z * b.m31,
					   a.x * b.m12 + a.y * b.m22 + a.z * b.m32,
					   a.x * b.m13 + a.y * b.m23 + a.z * b.m33};
	}

	vector3 operator*(const matrix3 &a, const vector3 &b)
	{
		return vector3(b.x * a.m11 + b.y * a.m12 + b.z * a.m13,
					   b.x * a.m21 + b.y * a.m22 + b.z * a.m23,
					   b.x * a.m31 + b.y * a.m32 + b.z * a.m33);
	}

	vector4 operator*(const vector4 &a, const matrix4 &b)
	{
		return vector4(a.x * b.m11 + a.y * b.m21 + a.z * b.m31 + a.w * b.m41,
					   a.x * b.m12 + a.y * b.m22 + a.z * b.m32 + a.w * b.m42,
					   a.x * b.m13 + a.y * b.m23 + a.z * b.m33 + a.w * b.m43,
					   a.x * b.m14 + a.y * b.m24 + a.z * b.m34 + a.w * b.m44);
	}

	vector4 operator*(const matrix4 &a, const vector4 &b)
	{
		return vector4(b.x * a.m11 + b.y * a.m12 + b.z * a.m13 + b.w * a.m14,
					   b.x * a.m21 + b.y * a.m22 + b.z * a.m23 + b.w * a.m24,
					   b.x * a.m31 + b.y * a.m32 + b.z * a.m33 + b.w * a.m34,
					   b.x * a.m41 + b.y * a.m42 + b.z * a.m43 + b.w * a.m44);
	}

	matrix4 RotationMatrixZ(float angle)
	{
		float cosb = cosf(angle);
		float sinb = sinf(angle);
		return matrix4{cosb, -sinb, 0, 0, sinb, cosb, 0, 0, 0, 0, 1.0f, 0, 0, 0, 0, 1.0f};
	}

	matrix4 RotationMatrixY(float angle)
	{
		float cosb = cosf(angle);
		float sinb = sinf(angle);
		return matrix4{cosb, 0, sinb, 0, 0, 1.0f, 0, 0, -sinb, 0, cosb, 0, 0, 0, 0, 1.0f};
	}

	matrix4 RotationMatrixX(float angle)
	{
		float cosb = cosf(angle);
		float sinb = sinf(angle);
		return matrix4{1.0f, 0, 0, 0, 0, cosb, -sinb, 0, 0, sinb, cosb, 0, 0, 0, 0, 1.0f};
	}

	matrix4 RotationMatrixAxis(const vector3 &origin, const vector3 &axis, float angle)
	{
		vector3 Ay = -axis.Normalize();

		vector3 Az;
		if (Ay == vector3(0, 0, 1.0f) || Ay == vector3(0, 0, -1.0f))
		{
			Az = (Ay ^ vector3(vector4(Ay, 1.0f) * RotationMatrixX(0.4f))).Normalize();
		}
		else
		{
			Az = (Ay ^ vector3(vector4(Ay, 1.0f) * RotationMatrixZ(0.4f))).Normalize();
		}
		vector3 Ax = Ay ^ Az;

		const vector3 &O = origin;

		float sinb = sinf(angle);
		float cosb = cosf(angle);

		float x1 = Ax.x * Ax.x * cosb + Ay.x * Ay.x + Az.x * Az.x * cosb;
		float y1 = Ax.y * (Ax.x * cosb + Az.x * sinb) + Ay.y * Ay.x + Az.y * (Az.x * cosb - Ax.x * sinb);
		float z1 = Ax.z * (Ax.x * cosb + Az.x * sinb) + Ay.z * Ay.x + Az.z * (Az.x * cosb - Ax.x * sinb);
		float w1 = O.x * (1 - Ax.x * Ax.x * cosb - Ay.x * Ay.x - Az.x * Az.x * cosb) +
				   O.y * (Az.y * (Ax.x * sinb - Az.x * cosb) + Ax.y * (-Ax.x * cosb - Az.x * sinb) - Ay.y * Ay.x) +
				   O.z * (Az.z * (Ax.x * sinb - Az.x * cosb) + Ax.z * (-Ax.x * cosb - Az.x * sinb) - Ay.z * Ay.x);

		float x2 = Ax.x * (Ax.y * cosb + Az.y * sinb) + Ay.x * Ay.y + Az.x * (Az.y * cosb - Ax.y * sinb);
		float y2 = Ax.y * Ax.y * cosb + Ay.y * Ay.y + Az.y * Az.y * cosb;
		float z2 = Ax.z * (Ax.y * cosb + Az.y * sinb) + Ay.z * Ay.y + Az.z * (Az.y * cosb - Ax.y * sinb);
		float w2 = O.y * (1 - Ax.y * Ax.y * cosb - Ay.y * Ay.y - Az.y * Az.y * cosb) +
				   O.x * (Az.x * (Ax.y * sinb - Az.y * cosb) + Ax.x * (-Ax.y * cosb - Az.y * sinb) - Ay.x * Ay.y) +
				   O.z * (Az.z * (Ax.y * sinb - Az.y * cosb) + Ax.z * (-Ax.y * cosb - Az.y * sinb) - Ay.z * Ay.y);

		float x3 = Ax.x * (Ax.z * cosb + Az.z * sinb) + Ay.x * Ay.z + Az.x * (Az.z * cosb - Ax.z * sinb);
		float y3 = Ax.y * (Ax.z * cosb + Az.z * sinb) + Ay.y * Ay.z + Az.y * (Az.z * cosb - Ax.z * sinb);
		float z3 = Ax.z * Ax.z * cosb + Ay.z * Ay.z + Az.z * Az.z * cosb;
		float w3 = O.z * (1 - Ax.z * Ax.z * cosb - Ay.z * Ay.z - Az.z * Az.z * cosb) +
				   O.x * (Az.x * (Ax.z * sinb - Az.z * cosb) + Ax.x * (-Ax.z * cosb - Az.z * sinb) - Ay.x * Ay.z) +
				   O.y * (Az.y * (Ax.z * sinb - Az.z * cosb) + Ax.y * (-Ax.z * cosb - Az.z * sinb) - Ay.y * Ay.z);

		return matrix4(x1, y1, z1, w1,
					   x2, y2, z2, w2,
					   x3, y3, z3, w3,
					   0, 0, 0, 1.0f);
	}

	matrix4 RotationMatrixRedirect(vector3 oldDir, vector3 newDir)
	{
		oldDir = oldDir.Normalize();
		newDir = newDir.Normalize();
		vector3 v = newDir ^ oldDir;
		if (v.Length() == 0)
			return IdentityMatrix;

		v = v.Normalize();
		vector3 o = oldDir ^ v;
		vector3 o2 = newDir ^ v;

		return matrix4(oldDir.x * newDir.x + o.x * o2.x + v.x * v.x, oldDir.y * newDir.x + o.y * o2.x + v.y * v.x, oldDir.z * newDir.x + o.z * o2.x + v.z * v.x, 0,
					   oldDir.x * newDir.y + o.x * o2.y + v.x * v.y, oldDir.y * newDir.y + o.y * o2.y + v.y * v.y, oldDir.z * newDir.y + o.z * o2.y + v.z * v.y, 0,
					   oldDir.x * newDir.z + o.x * o2.z + v.x * v.z, oldDir.y * newDir.z + o.y * o2.z + v.y * v.z, oldDir.z * newDir.z + o.z * o2.z + v.z * v.z, 0,
					   0, 0, 0, 1);
	}

	// Creates matrix to scale vector3 along x, y, z axes.
	matrix4 ScaleMatrix(float x, float y, float z)
	{
		return matrix4{x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1.0f};
	}

	matrix4 ScaleMatrixAxis(vector3 axis, float k)
	{
		axis = axis.Normalize();
		return matrix4(1 + axis.x * axis.x * (k - 1), axis.y * axis.x * (k - 1), axis.z * axis.x * (k - 1), 0,
					   axis.x * axis.y * (k - 1), 1 + axis.y * axis.y * (k - 1), axis.z * axis.y * (k - 1), 0,
					   axis.x * axis.z * (k - 1), axis.y * axis.z * (k - 1), 1 + axis.z * axis.z * (k - 1), 0,
					   0, 0, 0, 1);
	}

	// Creates matrix to offset vector3 along x, y, z axes.
	matrix4 OffsetMatrix(float offsetX, float offsetY, float offsetZ)
	{
		return matrix4{1.0f, 0, 0, offsetX, 0, 1.0f, 0.0f, offsetY, 0, 0, 1.0f, offsetZ, 0, 0, 0, 1.0f};
	}

	// Creates matrix to translate vector3 to camera's view space.
	matrix4 LookToMatrix(const vector3 &pos, const vector3 &dir, const vector3 &up)
	{
		vector3 Right = (up ^ dir).Normalize();
		vector3 Dir = dir.Normalize();
		vector3 Up = Dir ^ Right;

		return matrix4{Right.x, Right.y, Right.z, -pos.x * Right.x - pos.y * Right.y - pos.z * Right.z,
					   Up.x, Up.y, Up.z, -pos.x * Up.x - pos.y * Up.y - pos.z * Up.z,
					   Dir.x, Dir.y, Dir.z, -pos.x * Dir.x - pos.y * Dir.y - pos.z * Dir.z,
					   0, 0, 0, 1};
	}

	matrix4 LookAtMatrix(const vector3 &pos, const vector3 &at, const vector3 &up)
	{
		return LookToMatrix(pos, at - pos, up);
	}

	// Creates perspective projection matrix.
	matrix4 ProjectionMatrixPerspective(float fov, float nearPlane, float farPlane, float aspectRatio)
	{
		float tg = tanf(fov / 2);

		return matrix4{1.0f / (tg * aspectRatio), 0, 0, 0,
					   0, 1.0f / tg, 0, 0,
					   0, 0, farPlane / (farPlane - nearPlane), -nearPlane * farPlane / (farPlane - nearPlane),
					   0, 0, 1, 0};
	}

	matrix4 ProjectionMatrixOrthogonal(float width, float height, float farPlane)
	{
		return matrix4(2.0f / width, 0, 0, 0, 0, 2.0f / height, 0, 0, 0, 0, 1.0f / farPlane, 0, 0, 0, 0, 1.0f);
	}

	matrix4 ReflectMatrix(vector3 normal)
	{
		normal = normal.Normalize();
		return matrix4(1 - 2 * normal.x * normal.x, -2 * normal.y * normal.x, -2 * normal.z * normal.x, 0,
					   -2 * normal.x * normal.y, 1 - 2 * normal.y * normal.y, -2 * normal.z * normal.y, 0,
					   -2 * normal.x * normal.z, -2 * normal.y * normal.z, 1 - 2 * normal.z * normal.z, 0,
					   0, 0, 0, 1);
	}

	bool GenerateVertexTangents(const vector3 *vertices, const vector2 *texCoords, vector3 *tangents, vector3 *bitangents, int count)
	{
		if (count % 3 != 0)
			return false;
		for (int i = 0; i < count; i += 3)
		{
			matrix2 inversedUVs = matrix2(texCoords[i + 1].x - texCoords[i].x, texCoords[i + 1].y - texCoords[i].y,
										  texCoords[i + 2].x - texCoords[i].x, texCoords[i + 2].y - texCoords[i].y)
									  .Inverse();
			vector2 x = vector2(vertices[i + 1].x - vertices[i].x, vertices[i + 2].x - vertices[i].x) * inversedUVs;
			vector2 y = vector2(vertices[i + 1].y - vertices[i].y, vertices[i + 2].y - vertices[i].y) * inversedUVs;
			vector2 z = vector2(vertices[i + 1].z - vertices[i].z, vertices[i + 2].z - vertices[i].z) * inversedUVs;

			vector3 tang = vector3(x.x, y.x, z.x).Normalize();
			vector3 bitang = vector3(x.y, y.y, z.y).Normalize();
			tangents[i] = tang;
			tangents[i + 1] = tang;
			tangents[i + 2] = tang;
			bitangents[i] = bitang;
			bitangents[i + 1] = bitang;
			bitangents[i + 2] = bitang;
		}
		return true;
	}

	bool GenerateVertexTangentsIndexed(const vector3 *vertices, const vector2 *texCoords, vector3 *tangents, vector3 *bitangents, unsigned count, const unsigned *indices, unsigned indicesCount)
	{
		if (indicesCount % 3 != 0)
			return false;
		for (unsigned i = 0; i < indicesCount; i += 3)
		{
			matrix2 inversedUVs = matrix2(texCoords[indices[i + 1]].x - texCoords[indices[i]].x, texCoords[indices[i + 1]].y - texCoords[indices[i]].y,
										  texCoords[indices[i + 2]].x - texCoords[indices[i]].x, texCoords[indices[i + 2]].y - texCoords[indices[i]].y)
									  .Inverse();

			vector2 x = vector2(vertices[indices[i + 1]].x - vertices[indices[i]].x, vertices[indices[i + 2]].x - vertices[indices[i]].x) * inversedUVs;
			vector2 y = vector2(vertices[indices[i + 1]].y - vertices[indices[i]].y, vertices[indices[i + 2]].y - vertices[indices[i]].y) * inversedUVs;
			vector2 z = vector2(vertices[indices[i + 1]].z - vertices[indices[i]].z, vertices[indices[i + 2]].z - vertices[indices[i]].z) * inversedUVs;

			vector3 tang = vector3(x.x, y.x, z.x);
			vector3 bitang = vector3(x.y, y.y, z.y);
			tangents[indices[i]] += tang;
			tangents[indices[i + 1]] += tang;
			tangents[indices[i + 2]] += tang;
			bitangents[indices[i]] += bitang;
			bitangents[indices[i + 1]] += bitang;
			bitangents[indices[i + 2]] += bitang;
		}

		for (unsigned i = 0; i < count; i++)
		{
			tangents[i] = tangents[i].Normalize();
			bitangents[i] = bitangents[i].Normalize();
		}
		return true;
	}

	bool GenerateVertexNormalsCCW(const vector3 *vertices, vector3 *normals, unsigned count)
	{
		if (count % 3 != 0)
			return false;

		for (unsigned i = 0; i < count; i += 3)
		{
			normals[i] = ((vertices[i + 2] - vertices[i]) ^ (vertices[i + 1] - vertices[i])).Normalize();
			normals[i + 1] = normals[i];
			normals[i + 2] = normals[i];
		}

		return true;
	}

	bool GenerateVertexNormalsCW(const vector3 *vertices, vector3 *normals, unsigned count)
	{
		if (count % 3 != 0)
			return false;

		for (unsigned i = 0; i < count; i += 3)
		{
			normals[i] = ((vertices[i + 1] - vertices[i]) ^ (vertices[i + 2] - vertices[i])).Normalize();
			normals[i + 1] = normals[i];
			normals[i + 2] = normals[i];
		}

		return true;
	}

	bool GenerateVertexNormalsCCWIndexed(const vector3 *vertices, vector3 *normals, unsigned count, const unsigned *indices, unsigned indicesCount)
	{
		if (indicesCount % 3 != 0)
			return false;
		for (unsigned i = 0; i < indicesCount; i += 3)
		{
			vector3 normal = ((vertices[indices[i + 2]] - vertices[indices[i]]) ^ (vertices[indices[i + 1]] - vertices[indices[i]])).Normalize();
			normals[indices[i]] += normal;
			normals[indices[i + 1]] += normal;
			normals[indices[i + 2]] += normal;
		}

		for (unsigned i = 0; i < count; i++)
		{
			normals[i] = normals[i].Normalize();
		}

		return false;
	}

	bool GenerateVertexNormalsCWIndexed(const vector3 *vertices, vector3 *normals, unsigned count, const unsigned *indices, unsigned indicesCount)
	{
		if (indicesCount % 3 != 0)
			return false;
		for (unsigned i = 0; i < indicesCount; i += 3)
		{
			vector3 normal = ((vertices[indices[i + 1]] - vertices[indices[i]]) ^ (vertices[indices[i + 2]] - vertices[indices[i]])).Normalize();
			normals[indices[i]] += normal;
			normals[indices[i + 1]] += normal;
			normals[indices[i + 2]] += normal;
		}

		for (unsigned i = 0; i < count; i++)
		{
			normals[i] = normals[i].Normalize();
		}

		return false;
	}

	void GenerateSphericalTexCoords(const vector3 *vertices, vector2 *coords, unsigned count, float uOffset, float vOffset, bool uInverse, float uTile, float vTile)
	{
		for (unsigned i = 0; i < count; i++)
		{
			vector3 dir = vertices[i];
			dir.y /= dir.Length();
			float length = vector2(dir.x, dir.z).Length();
			dir.x /= length;
			dir.z /= length;

			coords[i].x = acosf(dir.x);
			if (dir.z < 0)
				coords[i].x = 2.0f * PI - coords[i].x;
			coords[i].x /= 2.0f * PI;
			if (uInverse)
				coords[i].x = 1 - coords[i].x;
			coords[i].y = acosf(dir.y) / PI;
			coords[i].x *= uTile;
			coords[i].y *= vTile;
			coords[i].x += uOffset;
			coords[i].y += vOffset;
		}
	}

} /* namespace math3d */
