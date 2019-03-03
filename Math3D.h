#pragma once

/*
NAMING CONVENTION

Static functions obey some naming conventions.

For matrices, function name consists of its meaning(essence), the keyword "Matrix" and the optional keyword.
Each part must start with an uppercase letter.
For example, "RotationMatrixAxis", where:
1) "Rotation" - the meaning of the function. It must rotate something;
2) "Matrix";
3) "Axis" - the optional keyword which says that we must define the rotation axis.
Another example, "ScaleMatrix", where:
1) "Scale" - means that we scale something, vector in this particular case.
2) "Matrix";
3) "" - it doesn't have the optional keyword.

For other functions this naming convention is the same except that we don't write "Matrix" keyword.

.....................................................................................................................
СОГЛАШЕНИЕ ОБ ИМЕНОВАНИИ

Статические функции подчиняются ряду правил об именах.

Для матриц имя функции состоит из ее значения (сути), ключевого слова «Matrix» и необязательного ключевого слова.
Каждая часть должна начинаться с буквы верхнего регистра.
Например, «RotationMatrixAxis», где:
1) «Rotation» - значение функции. Она должна что-то вращать;
2) «Matrix»;
3) «Axis» - необязательное ключевое слово, которое сообщает нам, что мы должны определить ось вращения.
Другой пример: «ScaleMatrix», где:
1) «Scale» - означает, что мы масштабируем что-то, вектор в этом конкретном случае.
2) «Matrix»;
3) "" - у имени «ScaleMatrix» нет необязательного ключевого слова.

Для других функций это соглашение об именах одинаково, за исключением того, что мы не пишем ключевое слово «Matrix».


*/

namespace Math3D
{
	const float PI = 3.1415926535897932384626433832795f;
	struct vector2;
	struct vector3;
	struct vector4;
	struct matrix2;
	struct matrix3;
	struct matrix4;

	struct vector2
	{
		union
		{
			struct
			{
				float x, y;
			};
			float all[2] = { 0 };
		};

		//Empty Constructor
		//Пустой конструктор
		vector2();
		//The constructor takes two parameters _x, _y for direct initialization
		//Конструктор принимает два параметра _x, _y для прямой инициализации
		vector2(float _x, float _y);
		//A constructor that takes a three-dimensional vector and initializes a two-dimensional vector with the first two coordinates
		//Конструктор который принимает трёхмерный вектор и инициализирует двумерный вектор двумя первыми координатами
		vector2(const vector3& a);
		//A constructor that takes a four-dimensional vector and initializes a two-dimensional vector with the first two coordinates
		//Конструктор который принимает четырёхмерный вектор и инициализирует двумерный вектор двумя первыми координатами
		vector2(const vector4& a);
		//The difference between the vectors A (Xa, Ya) and B (Xb, Yb). The result is a new vector C (Xa-Xb, Ya-Yb)
		//Разница векторов А(Xa,Ya) и B(Xb,Yb). Результатом будет новый вектор С(Xa-Xb,Ya-Yb)
		vector2 operator-(const vector2 & a) const;
		//Unary minus, all components of the vector A (Xa, Ya) change the sign to the opposite. The result is a new vector C (-Xa, -Ya)
		//Унарный минус, все компоненты вектора A(Xa,Ya) меняют знак на противоположный.Результатом будет новый вектор С(-Xa,-Ya)
		vector2 operator-() const;
		//Multiplication of the vector A (Xa, Ya) by the number C. The result is a new vector B (Xa * C, Ya * C)
		//Умножение вектора A(Xa,Ya) на число C.Результатом будет новый вектор B(Xa*C,Ya*C)
		vector2 operator*(float a) const;
		//Overloading the implicit reduction of a vector to an array in which the elements of the vector {x, y}
		//Перегрузка неявного приведения вектора к массиву в котором храняться элементы вектора {x,y}
		operator float*();
		//The sum of the product of the components of the vectors A (Xa, Ya) and B (Xb, Yb). The result is a float number Xa * Xb + Ya * Yb
		//Сумма произведения компонентов векторов А(Xa,Ya) и B(Xb,Yb). Результатом есть вещественное число Xa*Xb+Ya*Yb
		float operator*(const vector2 & a) const;
		//The division of the vector A (Xa, Ya) by the number C. The result is a new vector B (Xa / C, Ya / C)
		//Деление вектора А(Xa,Ya) на число С. Результатом будет новый вектор B(Xa/C,Ya/C)
		vector2 operator/(float a) const;
		//The length of the vector, which is calculated sqrtf (x * x + y * y)
		//Длина вектора, которая рассчитывается sqrtf(x*x + y*y)
		float Length() const;
		//The normalization of the vector A (Xa, Ya). Returns the normalized vector B (Xa / A.Length (), Ya / A.Length ())
		//Нормирование вектора А(Xa,Ya). Возвращает нормированный вектор B( Xa / А.Length() , Ya / А.Length() )
		vector2 Normalize() const;
		//The sum of the vector A (Xa, Ya) and the vector B (Xb, Yb). The result is a new vector C (Xa + Xb, Ya + Yb)
		//Сумма вектора А(Xa,Ya) и вектора B(Xb,Yb). Результатом будет новый вектор С(Xa+Xb,Ya+Yb)
		vector2 operator+(const vector2 & a) const;
		//The equality condition, if the vectors A (Xa, Ya) and B (Xb, Yb) have the SAME coordinates, then returns true
		//Условие равности, если вектора А(Xa,Ya) и B(Xb,Yb) имею ОДИНАКОВЫЕ кординаты, то возращает true
		bool operator==(const vector2 & a) const;
		//The sum of the vector A(Xa, Ya) and B(Xb, Yb). The result of the vector C(Xa + Xb, Ya + Yb) is written to the calling object
		//Сумма вектора А(Xa,Ya) и B(Xb,Yb). Результат вектор С(Xa+Xb,Ya+Yb) записывается в вызывающий объект
		vector2& operator+=(const vector2& a);
		//The difference between the vectors A (Xa, Ya) and B (Xb, Yb). The result of the vector C(Xa - Xb, Ya - Yb) is written to the calling object
		//Разница векторов А(Xa,Ya) и B(Xb,Yb). Результат вектор С(Xa-Xb,Ya-Yb) записывается в вызывающий объект
		vector2& operator-=(const vector2& a);
	};

	struct vector3
	{
		union
		{
			struct
			{
				float x, y, z;
			};
			float all[3] = { 0 };
		};
		//Empty Constructor
		//Пустой конструктор
		vector3();
		//The constructor takes three parameters _x, _y, _z for direct initialization
		//Конструктор принимает три параметра _x, _y, _z для прямой инициализации
		vector3(float _x, float _y, float _z);
		//A constructor that takes a two-dimensional vector and initializes a three-dimensional vector with the first two coordinates, and the third one is initialized with the additional parameter _z
		//Конструктор который принимает двумерный вектор и инициализирует трёхмерный вектор двумя первыми координатами, а третью кординату инициализирует дополнительным параметром _z
		vector3(const vector2& _xy, float _z);
		//A constructor that takes a four-dimensional vector and initializes a three-dimensional vector with the first three coordinates
		//Конструктор который принимает четырёхмерный вектор и инициализирует трёхмерный вектор тремя первыми координатами
		vector3(const vector4& _xyz);
		//The sum of the vector A (Xa, Ya, Za) and the vector B (Xb, Yb, Zb). The result is a new vector C (Xa + Xb, Ya + Yb, Za + Zb)
		//Сумма вектора А(Xa,Ya,Za) и вектора B(Xb,Yb,Zb). Результатом будет новый вектор С(Xa+Xb,Ya+Yb,Za+Zb)
		vector3 operator+(const vector3 & a) const;
		//The difference of the vector A (Xa, Ya, Za) and the vector B (Xb, Yb, Zb). The result is a new vector C (Xa-Xb, Ya-Yb, Za-Zb)
		//Разница вектора А(Xa,Ya,Za) и вектора B(Xb,Yb,Zb). Результатом будет новый вектор С(Xa-Xb,Ya-Yb,Za-Zb)
		vector3 operator-(const vector3 & a) const;
		//Unary minus, all components of the vector A (Xa, Ya, Za) change the sign to the opposite. The result is a new vector C (-Xa, -Ya, -Za)
		//Унарный минус, все компоненты вектора A(Xa,Ya,Za) меняют знак на противоположный.Результатом будет новый вектор С(-Xa,-Ya,-Za)
		vector3 operator-() const;
		//Multiplication of the vector A (Xa, Ya, Za) by the number C. The result is a new vector B (Xa * C, Ya * C, Za * C)
		//Умножение вектора A(Xa,Ya,Za) на число C. Результатом будет новый вектор B(Xa*C,Ya*C,Za*C)
		vector3 operator*(float a) const;
		//Overloading the implicit reduction of a vector to an array in which the elements of the vector {x, y, z}
		//Перегрузка неявного приведения вектора к массиву в котором храняться элементы вектора {x,y,z}
		operator float*();
		//The sum of the product of the components of the vectors A (Xa, Ya, Za) and B (Xb, Yb, Zb). The result is a float number Xa * Xb + Ya * Yb + Za * Zb
		//Сумма произведения компонентов векторов А(Xa,Ya,Za) и B(Xb,Yb,Zb). Результатом есть вещественное число Xa*Xb+Ya*Yb+Za*Zb
		float operator*(const vector3 & a) const;
		//The vector product of the vectors A (Xa, Ya, Za) and B (Xb, Yb, Zb). The result is a new vector C (Ya * Zb - Za * Yb, Za * Xb - Xa * Zb, Xa * Yb - Ya * Xb)
		//Векторное произведение векторов А(Xa,Ya,Za) и B(Xb,Yb,Zb). Результатом будет новый вектор С(Ya*Zb - Za*Yb, Za*Xb - Xa*Zb, Xa*Yb - Ya*Xb)
		vector3 operator^(const vector3 & a) const;
		//The division of the vector A (Xa, Ya, Za) by the number C. The result is a new vector B (Xa / C, Ya / C, Za / C)
		//Деление вектора А(Xa,Ya,Za) на число С. Результатом будет новый вектор B(Xa/C,Ya/C,Za/C)
		vector3 operator/(float a) const;
		//The length of the vector, which is calculated by sqrtf (x * x + y * y + z * z)
		//Длина вектора, которая расчитывается sqrtf(x*x + y*y + z*z)
		float Length() const;
		//The normalization of the vector A (Xa, Ya, Za). Returns the normalized vector B (Xa / A.Length (), Ya / A.Length (), Za / A.Length ())
		//Нормирование вектора А(Xa,Ya,Za). Возвращает нормированый вектор B( Xa / А.Length() , Ya / А.Length(),Za/А.Length() )
		vector3 Normalize() const;
		//The equality condition, if the vectors A (Xa,Ya,Za) and B (Xb,Yb,Zb) have the SAME coordinates, then returns true
		//Условие равности, если вектора А(Xa,Ya,Za) и B(Xb,Yb,Zb) имею ОДИНАКОВЫЕ кординаты, то возращает true
		bool operator==(const vector3 & a) const;
		//The sum of the vector A(Xa, Ya, Za) and B(Xb, Yb, Zb). The result of the vector C (Xa + Xb, Ya + Yb, Za+Zb) is written to the calling object
		//Сумма вектора А(Xa,Ya,Za) и B(Xb,Yb,Zb). Результат вектор С(Xa+Xb,Ya+Yb,Za+Zb) записывается в вызывающий объект
		vector3& operator+=(const vector3& a);
		//The difference between the vectors A(Xa,Ya,Za) and B(Xb,Yb,Zb). The result of the vector C(Xa-Xb,Ya-Yb,Za-Zb) is written to the calling object
		//Разница векторов А(Xa,Ya,Za) и B(Xb,Yb,Zb). Результат вектор С(Xa-Xb,Ya-Yb,Za-Zb) записывается в вызывающий объект
		vector3& operator-=(const vector3& a);
	};

	struct vector4
	{
		union
		{
			struct
			{
				float x, y, z, w;
			};
			float all[4] = { 0 };
		};
		//Empty Constructor
		//Пустой конструктор
		vector4();
		//The constructor takes the four parameters _x, _y, _z, _w for direct initialization
		//Конструктор принимает четырё параметра _x, _y, _z, _w для прямой инициализации
		vector4(float _x, float _y, float _z, float _w);
		//A constructor that takes a three-dimensional vector and initializes the four-dimensional vector with the first three coordinates, and the fourth one is initialized with the additional parameter _w
		//Конструктор который принимает трёхмерный вектор и инициализирует четырёхмерный вектор тремя первыми координатами, а четвёртую кординату инициализирует дополнительным параметром _w
		vector4(vector3 _xyz, float _w);
		//A constructor that takes a two-dimensional vector and initializes a four-dimensional vector with the first two coordinates, and the third and fourth coordinates are initialized with additional parameters _z and _w
		//Конструктор который принимает двумерный вектор и инициализирует четырёхмерный вектор двумя первыми координатами, а третью и четвёртую кординаты инициализирует дополнительными параметроми _z и _w
		vector4(vector2 _xy, float _z, float _w);
		//The sum of the vector A(Xa, Ya, Za, Wa) and the vector B(Xb, Yb, Zb, Wb). The result is a new vector C(Xa + Xb, Ya + Yb, Za + Zb, Wa + Wb)
		//Сумма вектора А(Xa,Ya,Za,Wa) и вектора B(Xb,Yb,Zb,Wb). Результатом будет новый вектор С(Xa+Xb,Ya+Yb,Za+Zb,Wa+Wb)
		vector4 operator+(const vector4 & a) const;
		//The difference between the vector A (Xa, Ya, Za, Wa) and the vector B(Xb, Yb, Zb, Wb). The result is a new vector C(Xa-Xb, Ya-Yb, Za-Zb, Wa-Wb)
		//Разница вектора А(Xa,Ya,Za,Wa) и вектора B(Xb,Yb,Zb,Wb). Результатом будет новый вектор С(Xa-Xb,Ya-Yb,Za-Zb,Wa-Wb)
		vector4 operator-(const vector4 & a) const;
		//Unary minus, all the components of the vector A(Xa, Ya, Za, Wa) change the sign to the opposite. The result is a new vector C(-Xa, -Ya, -Za, -Wa)
		//Унарный минус, все компоненты вектора A(Xa,Ya,Za,Wa) меняют знак на противоположный. Результатом будет новый вектор С(-Xa,-Ya,-Za,-Wa)
		vector4 operator-() const;
		//Multiplication of the vector A(Xa, Ya, Za, Wa) by the number C. The result is a new vector B(Xa * C, Ya * C, Za * C, Wa * C)
		//Умножение вектора A(Xa,Ya,Za,Wa) на число C. Результатом будет новый вектор B(Xa*C,Ya*C,Za*C,Wa*C)
		vector4 operator*(float a) const;
		//Overloading the implicit reduction of a vector to an array in which the elements of the vector {x, y, z, w}
		//Перегрузка неявного приведения вектора к массиву в котором храняться элементы вектора {x,y,z,w}
		operator float*();
		//The sum of the product of the components of the vectors A(Xa, Ya, Za, Wa) and B(Xb, Yb, Zb, Wb). The result is a float number Xa * Xb + Ya * Yb + Za * Zb + Wa * Wb
		//Сумма произведения компонентов векторов А(Xa,Ya,Za,Wa) и B(Xb,Yb,Zb,Wb). Результатом есть вещественное число Xa*Xb+Ya*Yb+Za*Zb+Wa*Wb
		float operator*(const vector4 & a) const;
		//The division of the vector A(Xa, Ya, Za, Wa) by the number C. The result is a new vector B(Xa / C, Ya / C, Za / C, Wa / C)
		//Деление вектора А(Xa,Ya,Za,Wa) на число С. Результатом будет новый вектор B(Xa/C,Ya/C,Za/C,Wa/C)
		vector4 operator/(float a) const;
		//The length of the vector, which is calculated by sqrtf (x * x + y * y + z * z + w * w)
		//Длина вектора, которая расчитывается sqrtf(x*x + y*y + z*z + w*w)
		float Length() const;
		//The normalization of the vector A(Xa, Ya, Za, Wa). Returns the normalized vector B(Xa / A.Length (), Ya / A.Length (), Za / A.Length (), Wa / A.Length ())
		//Нормирование вектора А(Xa,Ya,Za,Wa). Возвращает нормированый вектор B( Xa / А.Length() , Ya / А.Length(),Za/А.Length(),Wa/А.Length() )
		vector4 Normalize() const;
		//The equality condition, if the vectors А(Xa,Ya,Za,Wa) and B(Xb,Yb,Zb,Wb) have the SAME coordinates, then returns true
		//Условие равности, если вектора А(Xa,Ya,Za,Wa) и B(Xb,Yb,Zb,Wb) имею ОДИНАКОВЫЕ кординаты, то возращает true
		bool operator==(const vector4 & a) const;
		//The sum of the vector A(Xa,Ya,Za,Wa) and B(Xb,Yb,Zb,Wb). The result of the vector C(Xa+Xb,Ya+Yb,Za+Zb,Wa+Wb) is written to the calling object
		//Сумма вектора А(Xa,Ya,Za,Wa) и B(Xb,Yb,Zb,Wb). Результат вектор С(Xa+Xb,Ya+Yb,Za+Zb,Wa+Wb) записывается в вызывающий объект
		vector4& operator+=(const vector4& a);
		//The difference between the vectors A(Xa,Ya,Za,Wa) and B(Xb,Yb,Zb,Wb). The result of the vector C(Xa-Xb,Ya-Yb,Za-Zb,Wa-Wb) is written to the calling object
		//Разница векторов А(Xa,Ya,Za,Wa) и B(Xb,Yb,Zb,Wb). Результат вектор С(Xa-Xb,Ya-Yb,Za-Zb,Wa-Wb) записывается в вызывающий объект
		vector4& operator-=(const vector4& a);
	};

	struct matrix2
	{
		union
		{
			struct
			{
				float m11, m21, m12, m22;
			};
			float all[4] = { 0 };
			float all2D[2][2];
			vector2 columns[2];
		};
		//Empty Constructor
		//Пустой конструктор
		matrix2();
		//The constructor takes 4 parameters for direct matrix initialization
		//Конструктор принимает 4 параметра для прямой инициализации матрицы
		matrix2(float _m11, float _m21, float _m12, float _m22);
		//The constructor takes two two-dimensional vectors that initialize the matrix by columns
		//Конструктор принимает два двумерных вектора которые инициализируют матрицу по столбцам
		matrix2(const vector2& column1, const vector2& column2);
		//The constructor takes a 3x3 matrix and initializes its submatrix
		//Конструктор принимает матрицу 3x3 и инициализирует ее подматрицу
		matrix2(const matrix3& a, unsigned offsetX, unsigned offsetY);
		//The constructor takes a 4x4 matrix and initializes its submatrix
		//Конструктор принимает матрицу 4x4 и инициализирует ее подматрицу
		matrix2(const matrix4& a, unsigned offsetX, unsigned offsetY);
		//Finding the determinant of the matrix. Returns a float number
		//Нахождение детерминанта матрицы. Возвращает вещественное число
		float Determinant() const;
		//Finding the transpose matrix. Returns a new two-dimensional matrix
		//Нахождени транспонированной матрицы. Возвращает новую двумерную матрицу
		matrix2 Transpose() const;
		//Finding the inverse matrix. Returns a new two-dimensional matrix
		//Нахождение инверсной матрицы. Возвращает новую двумерную матрицу
		matrix2 Inverse() const;
		//Multiplication of a matrix by a number.Returns a new two - dimensional matrix
		//Умножение матрицы на число. Возвращает новую двумерную матрицу
		matrix2 operator*(float a) const;
		//Addition of two two-dimensional matrices. Returns a new two-dimensional matrix
		//Сложение двух двумерных матриц. Возвращает новую двумерную матрицу
		matrix2 operator+(const matrix2 & a) const;
		//The difference between two two-dimensional matrices. Returns a new two-dimensional matrix
		//Разница двух двумерных матриц. Возвращает новую двумерную матрицу
		matrix2 operator-(const matrix2 & a) const;
		//Multiplication of two two-dimensional matrices. Returns a new two-dimensional matrix
		//Умножение двух двумерных матриц. Возвращает новую двумерную матрицу
		matrix2 operator*(const matrix2 & a) const;
		//Access to a column by index
		//Доступ к столбцу по индексу
		vector2 & operator[](int index);

	};

	struct matrix3
	{
		union
		{
			struct
			{
				float m11, m21, m31, m12, m22, m32, m13, m23, m33;
			};
			float all[9] = { 0 };
			float all2D[3][3];
			vector3 columns[3];
		};
		//Empty Constructor
		//Пустой конструктор
		matrix3();
		//The constructor takes 9 parameters for direct matrix initialization
		//Конструктор принимает 9 параметров для прямой инициализации матрицы
		matrix3(float _m11, float _m21, float _m31, float _m12, float _m22, float _m32, float _m13, float _m23, float _m33);
		//The constructor takes three three-dimensional vectors that initialize the matrix by columns
		//Конструктор принимает три трехмерных вектора которые инициализируют матрицу по столбцам
		matrix3(const vector3& column1, const vector3& column2, const vector3& column3);
		//The constructor takes 2x2 matrix and initializes 3x3 matrix with 2x2 submatrix with xy offsets. Offsets mustn't be bigger than one!
		//Конструктор принимает матрицу 2x2 и инициализирует матрицу 3x3 подматрицей 2x2 с оступами по xy. Отступы должны быть не больше единицы!
		matrix3(const matrix2& a, unsigned offsetX, unsigned offsetY);
		//The constructor takes 4x4 matrix and initializes 3x3 matrix with 3x3 submatrix with xy offsets. Offsets mustn't be bigger than one!
		//Конструктор принимает матрицу 4x4 и инициализирует матрицу 3x3 подматрицей 3x3 с оступами по xy. Отступы должны быть не больше единицы!
		matrix3(const matrix4& a, unsigned offsetX, unsigned offsetY);
		//Finding the determinant of the matrix. Returns a float number
		//Нахождение детерминанта матрицы. Возвращает вещественное число
		float Determinant() const;
		//Finding the transpose matrix. Returns a new three-dimensional matrix
		//Нахождение транспонированной матрицы. Возвращает новую трёхмерную матрицу
		matrix3 Transpose() const;
		//Multiplication of a matrix by a number. Returns a new three-dimensional matrix
		//Умножение матрицы на число. Возвращает новую трёхмерную матрицу
		matrix3 operator*(float a) const;
		//Addition of two three-dimensional matrices. Returns a new three-dimensional matrix
		//Сложение двух трёхмерных матриц. Возвращает новую трёхмерную матрицу
		matrix3 operator+(const matrix3 & b) const;
		//The difference between two three-dimensional matrices. Returns a new three-dimensional matrix
		//Разница двух трёхмерных матриц. Возвращает новую трёхмерную матрицу
		matrix3 operator-(const matrix3 & b) const;
		//Multiplication of two three-dimensional matrices. Returns a new three-dimensional matrix
		//Умножение двух трёхмерных матриц. Возвращает новую трёхмерную матрицу
		matrix3 operator*(const matrix3 & b) const;
		//Finding the inverse matrix. Returns a new three-dimensional matrix
		//Нахождени инверсной матрицы. Возвращает новую трёхмерную матрицу
		matrix3 Inverse() const;
		//Overloading the implicit conversion of the matrix to an array of floats {m11, m21, m31,  m12, m22, m32,  m13, m23, m33}
		//Перегрузка неявного приведения матрицы к массиву в котором храняться элементы матрицы {m11, m21, m31,  m12, m22, m32,  m13, m23, m33}
		operator float*();
		//Access to a column by index
		//Доступ к столбцу по индексу
		vector3 & operator[](int index);
	};

	struct matrix4
	{
		union
		{
			struct
			{
				float m11, m21, m31, m41, m12, m22, m32, m42, m13, m23, m33, m43, m14, m24, m34, m44;
			};
			float all[16] = { 0 };
			float all2D[4][4];
			vector4 columns[4];
		};
		//Empty Constructor
		//Пустой конструктор
		matrix4();
		//The constructor takes 16 parameters for direct matrix initialization
		//Конструктор принимает 16 параметров для прямой инициализации матрицы
		matrix4(float _m11, float _m21, float _m31, float _m41, float _m12, float _m22, float _m32, float _m42, float _m13, float _m23, float _m33, float _m43, float _m14, float _m24, float _m34, float _m44);
		//The constructor takes 2x2 matrix and initializes 4x4 matrix with 2x2 submatrix with xy offsets. Offsets mustn't be bigger than two!
		//Конструктор принимает матрицу 2x2 и инициализирует матрицу 4x4 подматрицей 2x2 с оступами по xy. Отступы должны быть не больше двух!
		matrix4(const matrix2& a, unsigned offsetX, unsigned offsetY);
		//The constructor takes 3x3 matrix and initializes 4x4 matrix with 3x3 submatrix with xy offsets. Offsets mustn't be bigger than one!
		//Конструктор принимает матрицу 3x3 и инициализирует матрицу 4x4 подматрицей 3x3 с оступами по xy. Отступы должны быть не больше единицы!
		matrix4(const matrix3& a, unsigned offsetX, unsigned offsetY);
		//The constructor takes four four-dimensional vectors that initialize the matrix by columns
		//Конструктор принимает четырё четырёхмерных вектора которые инициализируют матрицу по столбцам
		matrix4(const vector4& column1, const vector4& column2, const vector4& column3, const vector4& column4);
		//Finding the determinant of the matrix. Returns a float number
		//Нахождение детерминанта матрицы. Возвращает вещественное число
		float Determinant() const;
		//Finding the transpose matrix. Returns a new four-dimensional matrix
		//Нахождение транспонированной матрицы. Возвращает новую четырёхмерную матрицу
		matrix4 Transpose() const;
		//Finding the inverse matrix. Returns a new four-dimensional matrix
		//Нахождение инверсной матрицы. Возвращает новую четырёхмерную матрицу
		matrix4 Inverse() const;
		//Multiplication of a matrix by a number. Returns a new four-dimensional matrix
		//Умножение матрицы на число.Возвращает новую четырёхмерную матрицу
		matrix4 operator*(float a) const;
		//Addition of two four-dimensional matrices. Returns a new four-dimensional matrix
		//Сложение двух четырёхмерных матриц. Возвращает новую четырёхмерную матрицу
		matrix4 operator+(const matrix4 & b) const;
		//The difference between two four-dimensional matrices. Returns a new four-dimensional matrix
		//Разница двух четырёхмерных матриц. Возвращает новую четырёхмерную матрицу
		matrix4 operator-(const matrix4 & b) const;
		//Multiplication of two four-dimensional matrices. Returns a new four-dimensional matrix
		//Умножение двух четырёхмерных матриц. Возвращает новую четырёхмерную матрицу
		matrix4 operator*(const matrix4 & b) const;
		//Overloading the implicit conversion of the matrix to an array of floats {m11, m21, m31, m41, m12, m22, m32, m42, m13, m23, m33, m43, m14, m24, m34, m44}
		//Перегрузка неявного приведения матрицы к массиву в котором храняться элементы матрицы {m11, m21, m31, m41, m12, m22, m32, m42, m13, m23, m33, m43, m14, m24, m34, m44}
		operator float*();
		//Access to a column by index
		//Доступ к столбцу по индексу
		vector4 & operator[](int index);
	};
	//Multiplication of a number by a four-dimensional matrix. Returns a new four-dimensional matrix
	//Умножение числа на четырёхмерную матрицу. Возвращает новую четырёхмерную матрицу
	matrix4 operator*(float a, const matrix4 & b);
	//Multiplication of a number by a three-dimensional matrix. Returns a new three-dimensional matrix
	//Умножение числа на трёхмерную матрицу. Возвращает новую трёхмерную матрицу
	matrix3 operator*(float a, const matrix3 & b);
	//Multiplication of a number by a two-dimensional matrix. Returns a new two-dimensional matrix
	//Умножение числа на двумерную матрицу. Возвращает новую двумерную матрицу
	matrix2 operator*(float a, const matrix2 & b);
	//Multiplication of a number by a three-dimensional vector. Returns a new three-dimensional vector
	//Умножение числа на трёхмерный вектор. Возвращает новый трёхмерный вектор
	vector3 operator*(float a, const vector3 & b);
	//Multiplication of a two-dimensional vector by a two-dimensional matrix. Returns a new two-dimensional vector
	//Умножение двумерного вектора на двумерную матрицу. Возвращает новый двумерный вектор
	vector2 operator*(const vector2 & a, const matrix2 & b);
	//Multiplication of a two-dimensional matrix by a two-dimensional vector. Returns a new two-dimensional vector
	//Умножение двумерной матрицы на двумерный вектор. Возвращает новый двумерный вектор
	vector2 operator*(const matrix2 & a, const vector2 & b);
	//Multiplication of a three-dimensional vector by a three-dimensional matrix. Returns a new three-dimensional vector
	//Умножение трёхмерного вектора на трёхмерную матрицу. Возвращает новый трёхмерный вектор
	vector3 operator*(const vector3 & a, const matrix3 & b);
	//Multiplication of a three-dimensional matrix by a three-dimensional vector. Returns a new three-dimensional vector
	//Умножение трёхмерной матрицы на трёхмерный вектор. Возвращает новый трёхмерный вектор
	vector3 operator*(const matrix3 & a, const vector3 & b);
	//Multiplication of a four-dimensional vector by a four-dimensional matrix. Returns a new four-dimensional vector
	//Умножение четырёхмерного вектора на четырёхмерную матрицу. Возвращает новый четырёхмерный вектор
	vector4 operator*(const vector4 & a, const matrix4 & b);
	//Multiplication of a four-dimensional matrix by a four-dimensional vector. Returns a new four-dimensional vector
	//Умножение четырёхмерной матрицы на четырёхмерный вектор. Возвращает новый четырёхмерный вектор
	vector4 operator*(const matrix4 & a, const vector4 & b);

	
	//Function creates rotation matrix to rotate 3D vector along the Z axis at an angle "angle".
	//If multiplying by a vector then the vector must stand leftward, be 4-dimensional and w coordinate must be equal to 1.
	//Функция создает матрицу вращения для поворота 3D вектора вдоль оси Z под углом «angle».
	//Если умножать на вектор, то вектор должен быть слева, быть 4-х мерным и w координата должна быть равна 1.
	matrix4 RotationMatrixZ(float angle);

	//Function creates rotation matrix to rotate 3D vector along the Y axis at an angle "angle".
	//If multiplying by a vector then the vector must stand leftward, be 4-dimensional and w coordinate must be equal to 1.
	//Функция создает матрицу вращения для поворота 3D вектора вдоль оси Y под углом «angle».
	//Если умножать на вектор, то вектор должен быть слева, быть 4-х мерным и w координата должна быть равна 1.
	matrix4 RotationMatrixY(float angle);

	//Function creates rotation matrix to rotate 3D vector along the X axis at an angle "angle".
	//If multiplying by a vector then the vector must stand leftward, be 4-dimensional and w coordinate must be equal to 1.
	//Функция создает матрицу вращения для поворота 3D вектора вдоль оси X под углом «angle».
	//Если умножать на вектор, то вектор должен быть слева, быть 4-х мерным и w координата должна быть равна 1.
	matrix4 RotationMatrixX(float angle);

	//Function creates rotation matrix to rotate 3D vector along the "axis" axis at an angle "angle".
	//If multiplying by a vector then the vector must stand leftward, be 4-dimensional and w coordinate must be equal to 1.
	//(Expansive operation!)
	//Функция создает матрицу вращения для поворота 3D вектора вдоль оси «axis» под углом «angle».
	//Если умножать на вектор, то вектор должен быть слева, быть 4-х мерным и w координата должна быть равна 1.
	//(Дорогостоящая операция!)
	matrix4 RotationMatrixAxis(const vector3& origin, const vector3& axis, float angle);

	//Function creates rotation matrix to redirect the 3D vector along the "newDir".
	//If multiplying by a vector then the vector must stand leftward, be 4-dimensional and w coordinate must be equal to 1.
	//Функция создаёт матрицу вращения для перенаправления 3D вектора вдоль "newDir".
	//Если умножать на вектор, то вектор должен быть слева, быть 4-х мерным и w координата должна быть равна 1.
	matrix4 RotationMatrixRedirect(vector3 oldDir, vector3 newDir);

	//Function creates matrix to scale 3D vector along x, y, z axes.
	//If multiplying by a vector then the vector must stand leftward, be 4-dimensional and w coordinate must be equal to 1.
	//Функция создает матрицу для масштабирования 3D вектора вдоль осей x, y, z.
	//Если умножать на вектор, то вектор должен быть слева, быть 4-х мерным и w координата должна быть равна 1.
	matrix4 ScaleMatrix(float x, float y, float z);

	//Function creates matrix to scale 3D vector along "Axis" axis.
	//If multiplying by a vector then the vector must stand leftward, be 4-dimensional and w coordinate must be equal to 1.
	//Функция создает матрицу для масштабирования трехмерного вектора вдоль оси «axis».
	//Если умножать на вектор, то вектор должен быть слева, быть 4-х мерным и w координата должна быть равна 1.
	matrix4 ScaleMatrixAxis(vector3 axis, float k);

	//Function creates matrix to offset 3D vector along x, y, z axes.
	//If multiplying by a vector then the vector must stand leftward, be 4-dimensional and w coordinate must be equal to 1.
	//Функция создает матрицу для смещения 3D-вектора вдоль осей x, y, z.
	//Если умножать на вектор, то вектор должен быть слева, быть 4-х мерным и w координата должна быть равна 1.
	matrix4 OffsetMatrix(float offsetX, float offsetY, float offsetZ);

	//Function creates matrix to translate 3D vector to camera's view space.
	//If multiplying by a vector then the vector must stand leftward, be 4-dimensional and w coordinate must be equal to 1.
	//Функция создает матрицу для преобразования 3D-вектора в пространство просмотра камеры.
	//Если умножать на вектор, то вектор должен быть слева, быть 4-х мерным и w координата должна быть равна 1.
	matrix4 LookToMatrix(const vector3 & pos, const vector3 & dir, const vector3 & up);

	//Function creates matrix to translate 3D vector to camera's view space.
	//If multiplying by a vector then the vector must stand leftward, be 4-dimensional and w coordinate must be equal to 1.
	//Функция создает матрицу для преобразования 3D-вектора в пространство просмотра камеры.
	//Если умножать на вектор, то вектор должен быть слева, быть 4-х мерным и w координата должна быть равна 1.
	matrix4 LookAtMatrix(const vector3 & pos, const vector3 & at, const vector3 & up);

	//Function creates perspective projection matrix. It projects vector on the "nearplane". 
	//If multiplying by a vector then the vector must stand leftward, be 4-dimensional and w coordinate must be equal to 1.
	//Функция создаёт матрицу для проектирования вектора на плоскость с учетом перспективы. Он проецирует вектор на «nearPlane».
	//Если умножать на вектор, то вектор должен быть слева, быть 4-х мерным и w координата должна быть равна 1.
	matrix4 ProjectionMatrixPerspective(float fov, float nearPlane, float farPlane, float aspectRatio);

	//Function creates orthogonal projection matrix.
	//If multiplying by a vector then the vector must stand leftward, be 4-dimensional and w coordinate must be equal to 1.
	//Функция создаёт матрицу для проектирования вектора на плоскость.
	//Если умножать на вектор, то вектор должен быть слева, быть 4-х мерным и w координата должна быть равна 1.
	matrix4 ProjectionMatrixOrthogonal(float width, float height, float farPlane);

	//Function creates matrix, multiplying by which reflects 3D vector from the plane defined by its "normal" vector.
	//If multiplying by a vector then the vector must stand leftward, be 4-dimensional and w coordinate must be equal to 1.
	//Функция создает матрицу, умножая на которую отражает 3D-вектор от плоскости, определяемой его «normal» вектором.
	//Если умножать на вектор, то вектор должен быть слева, быть 4-х мерным и w координата должна быть равна 1.
	matrix4 ReflectMatrix(vector3 normal);

	//Function generates tangent and bitangent vectors for each vertex in the array based on its texture coordinates. Vertex array mustn't be indexed!
	//Each next 3 vertices must form a new triangle so the total count of vertices must be multiple of 3.
	//Функция генерирует тангенциальные и бинагентные векторы для каждой вершины в массиве на основе ее координат текстуры. Массив Vertex не должен быть проиндексирован!
	//Каждые 3 вершины должны образовывать новый треугольник, поэтому общее количество вершин должно быть кратным 3.
	bool GenerateVertexTangents(const vector3* vertices, const vector2* texCoords, vector3* tangents, vector3* bitangents, int count);

	//Function generates tangent and bitangent vectors for each vertex in the array based on its texture coordinates and position. Vertex array must be indexed!
	//Each 3 indices must form a new triangle so the total count of indices must be multiple of 3.
	//Функция генерирует тангенциальные и бинагентные векторы для каждой вершины в массиве на основе ее координат и положения текстуры. Массив Vertex должен быть проиндексирован!
	//Каждые 3 индекса должны образовывать новый треугольник, поэтому общее количество индексов должно быть кратным 3.
	bool GenerateVertexTangentsIndexed(const vector3* vertices, const vector2* texCoords, vector3* tangents, vector3* bitangents, unsigned count, const unsigned* indices, unsigned indicesCount);

	//Function generates normal for each 3 vertices in the array based on their position.
	//Triangles must have the counterclockwise orientation to get the front facing normal.
	// Функция генерирует нормальное значение для каждой из трех вершин массива в зависимости от их положения.
	// Треугольники должны иметь ориентацию против часовой стрелки, чтобы получить вектор направленый вперед на нас.
	bool GenerateVertexNormalsCCW(const vector3* vertices, vector3* normals, unsigned count);

	//Function generates normal for each 3 vertices in the array based on their position.
	//Triangles must have the clockwise orientation to get the front facing normal.
	// Функция генерирует нормальное значение для каждой из трех вершин массива в зависимости от их положения.
	// Треугольники должны иметь ориентацию за часовой стрелкой, чтобы получить вектор направленый вперед на нас.
	bool GenerateVertexNormalsCW(const vector3* vertices, vector3* normals, unsigned count);

	//Function generates normal for each vertex in the array based on its position.
	//Each 3 indices must form a new triangle so the total count of indices must be multiple of 3.
	//Triangles must have the counterclockwise orientation to get the front facing normal.
	// Функция генерирует нормаль для каждой вершины в массиве на основе ее положения.
	// Каждые 3 индекса должны образовывать новый треугольник, поэтому общее количество индексов должно быть кратным 3.
	// Треугольники должны иметь ориентацию против часовой стрелки, чтобы получить вектор направленый вперед на нас.
	bool GenerateVertexNormalsCCWIndexed(const vector3* vertices, vector3* normals, unsigned count, const unsigned* indices, unsigned indicesCount);

	//Function generates normal for each vertex in the array based on its position.
	//Each 3 indices must form a new triangle so the total count of indices must be multiple of 3.
	//Triangles must have the clockwise orientation to get the front facing normal.
	// Функция генерирует нормаль для каждой вершины в массиве на основе ее положения.
	// Каждые 3 индекса должны образовывать новый треугольник, поэтому общее количество индексов должно быть кратным 3.
	// Треугольники должны иметь ориентацию за часовой стрелкой, чтобы получить вектор направленый вперед на нас.
	bool GenerateVertexNormalsCWIndexed(const vector3* vertices, vector3* normals, unsigned count, const unsigned* indices, unsigned indicesCount);

	//Function generates texture coordinates for mesh. Texture is projected spherically.
	//Функция генерирует координаты текстуры для меша. Текстура проецируется сферически.
	void GenerateSphericalTexCoords(const vector3* vertices, vector2* coords, unsigned count, float uOffset, float vOffset, bool uInverse, float uTile, float vTile);
	
	const matrix4 IdentityMatrix = matrix4(	1, 0, 0, 0,
											0, 1, 0, 0,
											0, 0, 1, 0,
											0, 0, 0, 1);
}