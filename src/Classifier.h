#pragma once

#include "Matrix.h"
#include "Vector.h"
#include "RenderingResources.h"

namespace SVM
{
	// TODO this shouldn't be in the SVM namespace right?
	class Classifier
	{
	public:
		virtual void Train(Matrix<double> const& TrainExamples, Vector<int> const& TrainLabels) = 0;
		virtual void TrainStep(Matrix<double> const& TrainExamples, Vector<int> const& TrainLabels, double step_size, int step_itters, double regularization) = 0;
		virtual Vector<int> Test(Matrix<double> const& TestExamples) = 0;
		virtual void Clear() = 0;
	};
}
