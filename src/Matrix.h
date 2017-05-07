#pragma once

#include "Vector.h"

template <class T>
class Matrix
{
public:
	Matrix();
	Matrix(int M, int N, bool columnMajor);
	Matrix(Matrix<T>&& oldMatrix);
	Matrix<T>& operator=(Matrix<T>&& oldMatrix);
	~Matrix();

	//Vector<T>& operator[](std::size_t i);

	static Vector<T> Multiply(Matrix<T> const& A, Vector<T> const& B);
	static Matrix<T> Multiply(Matrix<T> const& A, Matrix<T> const& B);
	static Matrix<T> Multiply(Matrix<T> const& A, T b);
	static Matrix<T> Add(Matrix<T> const& A, Matrix<T> const& B);
	static Matrix<T> Subtract(Matrix<T> const& A, Matrix<T> const& B);

	// Returns itself
	Matrix<T>& Multiply(T const& b);
	Matrix<T>& Add(Matrix<T> const& B);
	void AddAll(T const& b);

	Matrix<T> SubMatrix(int row_offset, int column_offset, int rows, int columns) const;
	Vector<T> GetRow(int i) const;
	Vector<T> GetColumn(int j) const;
	T GetElement(int i, int j) const;

	// T Determinant();
	// T Determinant(T sign, int i, int MM);

	void SetRow(int i, Vector<T> const& a);
	void SetColumn(int j, Vector<T> const& a);
	void SetElement(int i, int j, T x);

	void AddRow(int i, Vector<T> const& a);
	void AddColumn(int j, Vector<T> const& a);

	// Number of rows, length of a column
	int M;
	// Number of columns, length of a row
	int N;

	T* GetData() const;
private:
	Matrix(Matrix<T> const& copyMe) {};
	Matrix<T>& operator=(Matrix<T>& copyMe) {};

	bool columnMajor;
	T* data;
};
