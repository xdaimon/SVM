#include "pch.h"
#include "Matrix.h"
#include "Logger.h"
#define LOG(...)

using std::cout;
using std::endl;

template <class T>
Matrix<T>::Matrix() {
	data = nullptr;
	M = 0;
	N = 0;
}

template <class T>
Matrix<T>::Matrix(int M, int N, bool columnMajor) {
	LOG("%dx%d Matrix created", M, N);

	this->M = M;
	this->N = N;
	this->columnMajor = columnMajor;

	data = nullptr;

	// The matrix will remain uninitialized
	if (M * N == 0)
		return;

	data = (T*)malloc(M * N * sizeof(T));
	memset(data, 0, M * N * sizeof(T));
}

template <class T>
Matrix<T>::Matrix(Matrix<T>&& oldMatrix)
    : M(0),
      N(0),
      columnMajor(0),
      data(0) {
	LOG("Matrix move constructed");

	*this = std::move(oldMatrix);
}

template <class T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>&& oldMatrix) {
	LOG("Matrix move assigned");

	if (this == &oldMatrix)
		return oldMatrix;

	if (data != nullptr)
		free(data);

	M = oldMatrix.M;
	N = oldMatrix.N;
	data = oldMatrix.data;
	columnMajor = oldMatrix.columnMajor;

	oldMatrix.M = 0;
	oldMatrix.N = 0;
	oldMatrix.data = nullptr;
	oldMatrix.columnMajor = false;
	return *this;
}

template <class T>
Matrix<T>::~Matrix() {
	LOG("Matrix destroyed");

	if (data != nullptr)
		free(data);
}

template <class T>
Matrix<T> Matrix<T>::SubMatrix(int row_offset, int column_offset, int rows, int columns) const {
	if (row_offset + rows > M && column_offset + columns > N) {
		cout << "Out of bounds in SubMatrix" << endl;
		exit(-1);
	}
	Matrix<T> ret = Matrix<T>(rows, columns, columnMajor);
	for (int i = row_offset; i < row_offset + rows; ++i)
		for (int j = column_offset; j < column_offset + columns; ++j)
			ret.SetElement(i, j, GetElement(i, j));
	return ret;
}

template <class T>
Vector<T> Matrix<T>::Multiply(Matrix<T> const& A, Vector<T> const& B) {
	if (A.N != B.M) {
		LOG("Size of matrices submitted to Matrix::Multiply(Mat, Vec) invalid");
		return Vector<T>(0);
	}

	Vector<T> ret = Vector<T>(A.M);
	for (int i = 0; i < A.M; ++i)
		ret.SetComponent(i, Vector<T>::Dot(A.GetRow(i), B));
	return std::move(ret);
}

template <class T>
Matrix<T> Matrix<T>::Multiply(Matrix<T> const& A, Matrix<T> const& B) {
	if (A.N != B.M) {
		cout << "Size of matrices submitted to Matrix::Multiply(Mat, Mat) invalid" << endl;
		exit(-1);
		return Matrix<T>();
	}

	Matrix<T> ret = Matrix(A.M, B.N, false);
	for (int i = 0; i < A.M; ++i) {
		for (int k = 0; k < B.N; ++k) {
			double sum = 0.0;
			for (int j = 0; j < A.N; ++j)
				sum += A.GetElement(i, j) * B.GetElement(j, k);
			ret.SetElement(i, k, sum);
		}
	}
	return std::move(ret);
}

template <class T>
Matrix<T>& Matrix<T>::Multiply(T const& b) {
	for (int i = 0; i < M; ++i)
		for (int j = 0; j < N; ++j)
			SetElement(i, j, GetElement(i, j) * b);
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::Add(Matrix<T> const& B) {
	if (M != B.M || N != B.N) {
		cout << "Wrong sized matrix given to Matrix.Add(Mat)" << endl;
		exit(-1);
		return *this;
	}

	for (int i = 0; i < M; ++i)
		for (int j = 0; j < N; ++j)
			SetElement(i, j, GetElement(i, j) + B.GetElement(i, j));
	return *this;
}

template <class T>
Matrix<T> Matrix<T>::Multiply(Matrix<T> const& A, T b) {
	Matrix<T> ret = Matrix(A.M, A.N, A.columnMajor);
	for (int i = 0; i < A.M; ++i)
		for (int j = 0; j < A.N; ++j)
			ret.SetElement(i, j, A.GetElement(i, j) * b);
	return std::move(ret);
}

template <class T>
Matrix<T> Matrix<T>::Add(Matrix<T> const& A, Matrix<T> const& B) {
	if (A.M != B.M || A.N != B.N) {
		LOG("Incorrect matricies passed to Subtract");
		return Matrix();
	}
	Matrix<T> ret = Matrix(A.M, A.N, A.columnMajor);
	for (int i = 0; i < A.M; ++i)
		for (int j = 0; j < A.N; ++j)
			ret.SetElement(i, j, A.GetElement(i, j) + B.GetElement(i, j));
	return ret;
}

template <class T>
void Matrix<T>::AddAll(T const& b) {
	for (int i = 0; i < M; ++i)
		for (int j = 0; j < N; ++j)
			SetElement(i, j, GetElement(i, j) + b);
}

template <class T>
Matrix<T> Matrix<T>::Subtract(Matrix<T> const& A, Matrix<T> const& B) {
	if (A.M != B.M || A.N != B.N) {
		cout << "Incorrect matricies passed to Subtract" << endl;
		exit(-1);
		return Matrix();
	}
	Matrix<T> ret = Matrix(A.M, A.N, A.columnMajor);
	for (int i = 0; i < A.M; ++i)
		for (int j = 0; j < A.N; ++j)
			ret.SetElement(i, j, A.GetElement(i, j) - B.GetElement(i, j));
	return ret;
}

template <class T>
Vector<T> Matrix<T>::GetRow(int i) const {
	if (columnMajor) {
		cout << "GetRow not implemented for ColumnMajor mats" << endl;
		exit(-1);
	}

	if (data == nullptr) {
		LOG("Matrix data == nullptr");
		return Vector<T>(0);
	}

	if (i >= M) {
		LOG("Index out of bounds in GetRow");
		return Vector<T>(0);
	}

	return Vector<T>(N, data + i * N, false);
}

template <class T>
Vector<T> Matrix<T>::GetColumn(int j) const {
	if (!columnMajor) {
		cout << "GetColumn not implemneted on row major mats" << endl;
		;
		exit(-1);
	}

	if (data == nullptr) {
		LOG("Matrix data == nullptr");
		return Vector<T>(0);
	}

	if (j >= N) {
		LOG("Index out of bounds in GetRow");
		return Vector<T>(0);
	}

	return Vector<T>(M, data + M * j, false);
}

template <class T>
T Matrix<T>::GetElement(int i, int j) const {
	if (data == nullptr) {
		LOG("Matrix data == nullptr");
		return T(0);
	}

	if (i >= M || j >= N) {
		LOG("Index out of bounds in Matrix::GetElement");
		return T(0);
	}

	// if (i*j > M*N)
	// {
	// 	LOG("Index out of bounds in Matrix::GetElement");
	// 	return T(0);
	// }

	if (columnMajor)
		return data[i + j * M];
	else
		return data[i * N + j];
}

template <class T>
void Matrix<T>::SetRow(int i, Vector<T> const& a) {
	if (data == nullptr)
		return LOG("Matrix data == nullptr");

	if (columnMajor)
		return LOG("SetRow to columnMajor matrix unimplemented");

	if (i >= M)
		return LOG("Index out of bounds in SetRow");

	memcpy(data + i * N, a.GetData(), N * sizeof(T));
}

template <class T>
void Matrix<T>::SetColumn(int j, Vector<T> const& a) {
	if (!columnMajor)
		return LOG("SetColumn to non columnMajor matrix unimplemented");

	if (j >= N)
		return LOG("Index out of bounds in SetColumn");

	if (data == nullptr)
		return LOG("Matrix data == nullptr");

	memcpy(data + M * j, a.GetData(), M * sizeof(T));
}

template <class T>
void Matrix<T>::SetElement(int i, int j, T x) {
	if (data == nullptr)
		return LOG("Matrix data == nullptr");

	if (columnMajor)
		data[i + j * M] = x;
	else
		data[i * N + j] = x;
}

template <class T>
void Matrix<T>::AddRow(int i, Vector<T> const& a) {
	if (columnMajor)
		return LOG("AddRow to columnMajor matrix unimplemented");

	if (i > M)
		return LOG("Index out of bounds in AddRow");

	// Not yet initialized
	if (M * N == 0)
		N = a.M;

	if (M * N != 0)
		if (a.M != N)
			return LOG("Wrong size vector given to Matrix::AddRow");

	data = (T*)realloc(data, (M + 1) * N * sizeof(T));
	memmove(data + (i + 1) * N, data + i * N, (M - i) * N * sizeof(T));
	M++;
	SetRow(i, a);
}

template <class T>
void Matrix<T>::AddColumn(int j, Vector<T> const& a) {
	if (!columnMajor)
		return LOG("AddColumn to non columnMajor matrix unimplemented");

	if (j > N)
		return LOG("Index out of bounds in AddColumn");

	// Not yet initialized
	if (M * N == 0)
		M = a.M;

	if (M * N != 0)
		if (a.M != M)
			return LOG("Wrong size vector given to Matrix::AddColumn");

	data = (T*)realloc(data, M * (N + 1) * sizeof(T));
	memmove(data + M * (j + 1), data + M * j, M * (N - j) * sizeof(T));
	N++;
	SetColumn(j, a);
}

template <class T>
T* Matrix<T>::GetData() const {
	if (data == nullptr)
		LOG("Matrix data == nullptr");
	return data;
}

template class Matrix<double>;
template class Matrix<int>;
