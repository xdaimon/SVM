#include "pch.h"
#include "RenderingResources.h"
#include "Logger.h"

SVM::RenderingResources::RenderingResources() {
	LOG("RenderingResources created");
}

SVM::RenderingResources::RenderingResources(RenderingResources&& other) {
	LOG("RenderingResources move created");

	data = std::move(other.data);
	fragmentShader = std::move(other.fragmentShader);
}

SVM::RenderingResources& SVM::RenderingResources::operator=(RenderingResources&& other) {
	LOG("RenderingResources move assigned");
	if (this == &other)
		return other;

	data = std::move(other.data);
	fragmentShader = std::move(other.fragmentShader);
	return *this;
}

SVM::RenderingResources::~RenderingResources() {
	LOG("RenderingResources destructed");
}

SVM::ResourceArray::ResourceArray()
    : data(nullptr) {
	LOG("ResourceArray created");
}

SVM::ResourceArray::ResourceArray(ResourceArray&& other) {
	LOG("ResourceArray move created");
	data = other.data;
	ownsData = other.ownsData;
	initializer = std::move(other.initializer);
	uploader = std::move(other.uploader);

	other.data = nullptr;
}

SVM::ResourceArray& SVM::ResourceArray::operator=(ResourceArray&& other) {
	LOG("ResourceArray move assigned");
	if (this == &other)
		return other;

	if (ownsData)
		if (data != nullptr)
			free(data);

	data = other.data;
	ownsData = other.ownsData;
	initializer = std::move(other.initializer);
	uploader = std::move(other.uploader);

	other.data = nullptr;

	return *this;
}

SVM::ResourceArray::~ResourceArray() {
	LOG("RenderingResources destructed");
	if (ownsData)
		if (data != nullptr)
			free(data);
}
