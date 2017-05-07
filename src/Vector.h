#pragma once

template <class T>
class Vector
{
public:
	Vector();
	Vector(int M);
	Vector(int M, T* data, bool ownsMemory);
	Vector(Vector<T>&& oldVector);
	Vector<T>& operator=(Vector<T>&& oldVector);
	~Vector();

	//T& operator[](int i);

	static T Dot(Vector<T> const& a, Vector<T> const& b);
	static Vector<T> Multiply(Vector<T> const& a, T b);
	static Vector<T> Sum(Vector<T> const& a, Vector<T> const& b);
	static Vector<T> Subtract(Vector<T> const& a, Vector<T> const& b);
	static T Length(Vector<T> const& a);

	Vector<T> SubVector(int start, int components) const;

	T GetComponent(int i) const;

	void SetComponent(int i, T x);
	void Append(T x);

	void AddComponent(int i, T x);

	int M;

	T* GetData() const;
private:
	Vector(const Vector<T>& copyMe) {};
	Vector<T>& operator=(Vector<T>& copyMe) { };

	bool ownsMemory;
	T* data;
};
