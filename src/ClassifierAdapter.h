#pragma once

#include "pch.h"

#include "MulticlassSVM.h"
#include "Matrix.h"
#include "Vector.h"

namespace SVM
{
	enum class RenderState {
		Render3D,
		Render2D
	};

	class ClassifierAdapter
	{
	public:
		virtual void Setup(SVM::RenderState state) = 0;
		virtual void Update() = 0;
		virtual void Train() = 0;
		virtual void TrainStep(double step_size, int step_itters, double regularization) = 0;
		virtual void PrintInformation() = 0;
		// virtual void WriteInformation() = 0;
		virtual void Clear() = 0;

		virtual void AddTrainExample(Vector<double> const& example, int const label) = 0;
		virtual void AddTestExample(Vector<double> const& example, int const label) = 0;
		virtual void OnInputChange(Vector<double> positionInfo) = 0;

		virtual void SetDataSet(int i) = 0;

		virtual RenderingResources& GetRenderingResources() = 0;
		virtual Vector<double> GetClassColor(int i) = 0;
	};
}
