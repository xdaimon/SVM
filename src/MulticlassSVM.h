#pragma once

#include "Classifier.h"
#include "Matrix.h"
#include "Vector.h"
#include "RenderingResources.h"

namespace SVM
{
	static const int MaxNumClasses = 6;
	static const int MaxNumExamples = 300;
	static const int NumDataSets = 2;

	class MulticlassSVM : public Classifier
	{
	public:
		// Inherited from Classifier interface
		void Train(Matrix<double> const& TrainExamples, Vector<int> const& TrainLabels);
		void TrainStep(Matrix<double> const& TrainExamples, Vector<int> const& TrainLabels, double step_size, int step_itters, double regularization);
		Vector<int> Test(Matrix<double> const& TestExampless);
		void Clear();

		Matrix<double> Weights;
		Vector<int> RowToClassMap;
	private:
		double weights_regularization = 1.;
		void FindLabelsInUse(Vector<int> const& TrainLabels);
		int trainStepItter;
		double Loss(Matrix<double> const & trainExamples, Vector<int> const & trainLabels, Vector<int> const & classToRowMap);
		double previous_loss = 1.;
		Matrix<double> LossGradient(Matrix<double> const & trainExamples, Vector<int> const & trainLabels, Vector<int> const & classToRowMap);
	};
}
