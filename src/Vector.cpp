#include "pch.h"
#include "Vector.h"
#include "Logger.h"
#define LOG(...)

template <class T>
Vector<T>::Vector() {
	data = nullptr;
	M = 0;
}

template <class T>
Vector<T>::Vector(int M) {
	this->M = M;

	data = nullptr;

	if (M != 0) {
		data = (T*)malloc(M * sizeof(T));
		memset(data, 0, M * sizeof(T));
	}

	LOG("%dD Vector created which owns data at %x", M, data);

	ownsMemory = true;
}

template <class T>
Vector<T>::Vector(Vector<T>&& oldVector)
    : M(0),
      ownsMemory(0),
      data(0) {
	LOG("Vector move constructed");

	*this = std::move(oldVector);
}

template <class T>
Vector<T>& Vector<T>::operator=(Vector<T>&& oldVector) {
	LOG("Vector move assigned");

	if (this == &oldVector)
		return oldVector;

	if (ownsMemory && data != nullptr)
		free(data);

	M = oldVector.M;
	data = oldVector.data;
	ownsMemory = oldVector.ownsMemory;

	oldVector.M = 0;
	oldVector.data = nullptr;
	oldVector.ownsMemory = false;

	return *this;
}

template <class T>
Vector<T>::Vector(int M, T* data, bool ownsMemory) {
	if (ownsMemory)
		LOG("%dD Vector created which owns data at %x", M, data);
	else
		LOG("%dD Vector created which shares data %x", M, data);

	this->M = M;
	this->data = data;
	this->ownsMemory = ownsMemory;
}

template <class T>
Vector<T>::~Vector() {
	LOG("Vector destroyed");

	if (ownsMemory)
		if (data != nullptr)
			free(data);
}

/*
template<class T>
T & Vector<T>::operator[](int i) {
	if (i > M) {
		LOG("Index out of bounds in Vector::GetComponent");
		return *data;
	}
	if (data == nullptr) {
		LOG("Vector data == nullptr");
		exit(-1);
		return *data;
	}

	return *(data+i);
}
*/

template <class T>
T Vector<T>::Dot(Vector<T> const& a, Vector<T> const& b) {
	if (a.M != b.M) {
		LOG("Size of vectors submitted to Vector::Dot do not match");
		exit(-1);
		return T(0);
	}

	T sum = T(0);
	for (int i = 0; i < a.M; ++i)
		sum += a.GetComponent(i) * b.GetComponent(i);
	LOG("asdf %f ", sum);
	return sum;
}

template <class T>
Vector<T> Vector<T>::Multiply(Vector<T> const& a, T b) {
	Vector<T> ret = Vector<T>(a.M);
	for (int i = 0; i < a.M; ++i)
		ret.SetComponent(i, a.GetComponent(i) * b);
	return std::move(ret);
}

template <class T>
Vector<T> Vector<T>::Sum(Vector<T> const& a, Vector<T> const& b) {
	Vector<T> ret = Vector<T>(a.M);
	for (int i = 0; i < a.M; ++i)
		ret.SetComponent(i, a.GetComponent(i) + b.GetComponent(i));
	return ret;
}

template <class T>
Vector<T> Vector<T>::Subtract(Vector<T> const& a, Vector<T> const& b) {
	Vector<T> ret = Vector<T>(a.M);
	for (int i = 0; i < a.M; ++i)
		ret.SetComponent(i, a.GetComponent(i) - b.GetComponent(i));
	return ret;
}

template <class T>
T Vector<T>::Length(Vector<T> const& a) {
	T ret = T(0);
	for (int i = 0; i < a.M; ++i)
		ret += a.GetComponent(i) * a.GetComponent(i);
	return ret;
}

template <class T>
Vector<T> Vector<T>::SubVector(int start, int components) const {
	if (start + components > M) {
		LOG("Out of bounds in SubVector");
		return Vector<T>();
	}
	Vector<T> ret = Vector<T>(components);
	for (int i = start; i < start + components; ++i)
		ret.SetComponent(i, GetComponent(i));
	return ret;
}

template <class T>
void Vector<T>::SetComponent(int i, T x) {
	if (i >= M)
		return LOG("Index out of bounds in Vector::SetComponent");
	if (data == nullptr)
		return LOG("Vector data == nullptr");

	data[i] = x;
}

template <class T>
void Vector<T>::Append(T x) {
	data = (T*)realloc(data, (M + 1) * sizeof(T));
	data[M] = x;
	M++;
}

template <class T>
void Vector<T>::AddComponent(int i, T x) {
	if (i > M)
		return LOG("Index out of bounds in Vector::AddComponent");

	data = (T*)realloc(data, (M + 1) * sizeof(T));
	memmove(data + i + 1, data + i, (M - i) * sizeof(T));
	M++;
	data[i] = x;
}

template <class T>
T Vector<T>::GetComponent(int i) const {
	if (i > M) {
		LOG("Index out of bounds in Vector::GetComponent");
		return T(0);
	}
	if (data == nullptr) {
		LOG("Vector data == nullptr");
		return T(0);
	}

	return data[i];
}

template <class T>
T* Vector<T>::GetData() const {
	if (data == nullptr)
		LOG("Vector data == nullptr");
	return data;
}

template class Vector<double>;
template class Vector<int>;
