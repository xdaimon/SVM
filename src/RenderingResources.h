#pragma once

#include <vector>
#include <string>
#include <functional>

namespace SVM
{
	class ResourceArray
	{
	public:
		ResourceArray();
		ResourceArray(ResourceArray&& oldMatrix);
		ResourceArray& operator=(ResourceArray&& other);
		~ResourceArray();

		float* data;
		bool ownsData;

		std::function<void(ResourceArray & ra)> initializer;
		std::function<void(ResourceArray const& ra)> uploader;

	private:
		ResourceArray(ResourceArray const& copyMe) {};
		ResourceArray& operator=(ResourceArray& copyMe) {};
	};

	class RenderingResources
	{
	public:
		RenderingResources();
		RenderingResources(RenderingResources&& oldMatrix);
		RenderingResources& operator=(RenderingResources&& other);
		~RenderingResources();

		std::string fragmentShader;
		std::vector<ResourceArray> data;

	private:
		RenderingResources(RenderingResources const& copyMe) {};
		RenderingResources& operator=(RenderingResources& copyMe) {};
	};
}
