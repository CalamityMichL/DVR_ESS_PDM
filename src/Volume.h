/* Copyright Michael Rauter
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * Licensed under the Apache License, Version 2.0 the "License";
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include "SimpleITK.h"
#include "sitkImageFileReader.h"
#include "sitkGreaterEqualImageFilter.h"
#include "sitkLessEqualImageFilter.h"
#include "sitkAndImageFilter.h"
#include "sitkStatisticsImageFilter.h"

#include <iomanip>

namespace sitk = itk::simple;

class Volume {
private:
	class ProgressUpdate
		: public sitk::Command
	{
	public:
		ProgressUpdate(const sitk::ProcessObject& po, const Volume& refV)
			: m_Process(po), refVolume(refV)
		{}

		void Execute() override
		{
			refVolume.PropagateLoadProgress(m_Process.GetName(), m_Process.GetProgress());
		}
	private:
		const sitk::ProcessObject& m_Process;
		const Volume& refVolume;
	};

	void (*callbackOnLoadProgress)(const std::string caption, float progress) = nullptr;
	void (*callbackLoadProgressStart)() = nullptr;
	void (*callbackLoadProgressStop)() = nullptr;

public:
	sitk::Image *volume;	
	sitk::Image *volumeFloat;
	sitk::Image volumeNative;
	sitk::Image volumeConverted;
	sitk::Image volumeUint16;

	double originalMinValue = 0.0;
	double originalMaxValue = 0.0;

	double minIntensityValue = 0.0;
	double maxIntensityValue = 0.0;

	Volume(	const std::vector<std::string>& inputImageFileNames,
			const std::string& strImageIOType = "",
			void (*funCallbackOnLoadProgress)(const std::string caption, float progress) = nullptr,
			void (*funCallbackLoadProgressStart)() = nullptr,
			void (*funCallbackLoadProgressStop)() = nullptr
		)
	{
		if (funCallbackLoadProgressStart != nullptr)
		{
			callbackLoadProgressStart = funCallbackLoadProgressStart;
			funCallbackLoadProgressStart();
		}
		if (funCallbackLoadProgressStop != nullptr)
		{
			callbackLoadProgressStop = funCallbackLoadProgressStop;
		}
		if (funCallbackOnLoadProgress != nullptr)
		{
			callbackOnLoadProgress = funCallbackOnLoadProgress;
		}

		std::unique_ptr<sitk::ImageReaderBase> reader;
		std::unique_ptr <ProgressUpdate> progressUpdate;

		if (inputImageFileNames.size() == 1)
		{
			reader = std::make_unique<sitk::ImageFileReader>();
			sitk::ImageFileReader* imageFileReader = dynamic_cast<sitk::ImageFileReader*>(reader.get());
			imageFileReader->SetFileName(inputImageFileNames[0]);
			imageFileReader->SetOutputPixelType(sitk::sitkUnknown);
			imageFileReader->SetImageIO("");
			sitk::ProcessObject *po = dynamic_cast<sitk::ProcessObject*>(imageFileReader);
			progressUpdate = std::make_unique<ProgressUpdate>(*po, *this);
			imageFileReader->AddCommand(sitk::sitkProgressEvent, *progressUpdate);
		}
		else
		{
			reader = std::make_unique<sitk::ImageSeriesReader>();
			sitk::ImageSeriesReader* imageSeriesReader = dynamic_cast<sitk::ImageSeriesReader*>(reader.get());
			imageSeriesReader->SetFileNames(inputImageFileNames);
			sitk::ProcessObject* po = dynamic_cast<sitk::ProcessObject*>(imageSeriesReader);
			progressUpdate = std::make_unique<ProgressUpdate>(*po, *this);
			imageSeriesReader->AddCommand(sitk::sitkProgressEvent, *progressUpdate);
		}
		volumeNative = reader->Execute();

		std::vector<unsigned int> padding = { 0,0,0 };
		std::vector<unsigned int> initialVolumeDims = volumeNative.GetSize();
		for (unsigned int i = 0; i < volumeNative.GetDimension(); i++)
			padding[i] = initialVolumeDims[i] % 2;

		if (padding[0] > 0 || padding[1] > 0 || padding[2] > 0)
			volumeNative = sitk::ConstantPad(volumeNative, { 0,0,0 }, padding, 0.0);

		sitk::MinimumMaximumImageFilter minMaxFilter = sitk::MinimumMaximumImageFilter();
		minMaxFilter.Execute(volumeNative);
		double highestIntensity = minMaxFilter.GetMaximum();
		double lowestIntensity = minMaxFilter.GetMinimum();

		volume = &volumeNative;
		std::vector<std::string> metadata = volume->GetMetaDataKeys();
		originalMinValue = lowestIntensity;
		originalMaxValue = highestIntensity;
		minIntensityValue = lowestIntensity;
		maxIntensityValue = highestIntensity;

		itk::simple::PixelIDValueEnum pixelID = GetPixelID();
		if (pixelID == itk::simple::PixelIDValueEnum::sitkFloat32 || pixelID == itk::simple::PixelIDValueEnum::sitkFloat64)
		{
			volumeFloat = &volumeNative;

			sitk::SubtractImageFilter subtractImageFilter = sitk::SubtractImageFilter();
			*volumeFloat = subtractImageFilter.Execute(*volumeFloat, lowestIntensity);

			double intensityRangeSpan = highestIntensity - lowestIntensity + 1;

			double factor = (double)(UINT16_MAX) / intensityRangeSpan;
			sitk::MultiplyImageFilter multiplyFilter = sitk::MultiplyImageFilter();
			*volumeFloat = multiplyFilter.Execute(*volumeFloat, factor);

			volumeUint16 = sitk::Cast(*volumeFloat, sitk::PixelIDValueEnum::sitkUInt16);

			volume = &volumeUint16;

			minIntensityValue = 0.0;
			maxIntensityValue = (double)(UINT16_MAX);
		}

		if (funCallbackLoadProgressStop != nullptr)
		{
			funCallbackLoadProgressStop();
		}
	}

	void PermuteAxes(std::vector<unsigned int> order) {
		*volume = sitk::PermuteAxes(*volume, order);
	}

	void Flip(std::vector<bool> flipAxes, bool flipAboutOrigin = false) {
		*volume = sitk::Flip(*volume, flipAxes, flipAboutOrigin);
	}

	void GetHistogram2(std::vector<uint32_t>& hist, std::vector<uint32_t>& edges)
	{
		hist.clear();
		edges.clear();

		int numBins = 256;

		edges.resize(numBins + 1);
		for (int i = 0; i < numBins + 1; i++)
			edges[i] = i;

		hist.resize(numBins);

		for (int i = 0; i < 256; i++)
		{
			sitk::GreaterEqualImageFilter ge_filter = sitk::GreaterEqualImageFilter();
			sitk::LessEqualImageFilter le_filter = sitk::LessEqualImageFilter();
			double lower = (minIntensityValue + (i/256.0) * (maxIntensityValue - minIntensityValue + 1));
			double upper = (minIntensityValue + ((i+1)/256.0) * (maxIntensityValue - minIntensityValue + 1));
			sitk::Image res = ge_filter.Execute(*volume, lower);
			sitk::Image res2 = le_filter.Execute(*volume, upper);

			sitk::AndImageFilter and_filter = sitk::AndImageFilter();
			res = and_filter.Execute(res, res2);


			sitk::StatisticsImageFilter statisticsFilter = sitk::StatisticsImageFilter();
			statisticsFilter.Execute(res);
			double histCount = statisticsFilter.GetSum();

			hist[i] = (uint32_t)histCount;
		}
	}

	void GetHistogram(std::vector<uint32_t> &hist, std::vector<uint32_t> &edges)
	{
		hist.clear();
		edges.clear();

		int numBins = 256;
		
		edges.resize(numBins + 1);
		for (int i = 0; i < numBins + 1; i++)
			edges[i] = i;

		hist.resize(numBins);

		if (callbackLoadProgressStart != nullptr)
			callbackLoadProgressStart();		

		int progress_percentage = 0;
		callbackOnLoadProgress("computing voxel histogram", progress_percentage);
		for (int i = 0; i < volume->GetNumberOfPixels(); i++)
		{
			int current_progress = (int)((float)i / (float)volume->GetNumberOfPixels() * 100.0f);
			if (current_progress > progress_percentage)
			{
				progress_percentage = current_progress;
				callbackOnLoadProgress("computing voxel histogram", progress_percentage);
			}
			
			uint32_t idx;
			double value;
			switch (GetPixelID()) {
			case itk::simple::PixelIDValueEnum::sitkUInt8:
				value = (double)(volume->GetBufferAsUInt8()[i]);
				idx = (uint32_t)(value);
				break;
			case itk::simple::PixelIDValueEnum::sitkInt8:
				value = (double)(volume->GetBufferAsInt8()[i]);
				idx = (uint32_t)((value - minIntensityValue) / (maxIntensityValue - minIntensityValue + 1) * 256);
				break;
			case itk::simple::PixelIDValueEnum::sitkUInt16:
				value = (double)(volume->GetBufferAsUInt16()[i]);
				idx = (uint32_t)((value - minIntensityValue) / (maxIntensityValue - minIntensityValue + 1) * 256);
				break;
			case itk::simple::PixelIDValueEnum::sitkInt16:
				value = (double)(volume->GetBufferAsInt16()[i]);
				idx = (uint32_t)((value - minIntensityValue) / (maxIntensityValue - minIntensityValue + 1) * 256);
				break;
			case itk::simple::PixelIDValueEnum::sitkFloat32:
				value = (double)(volume->GetBufferAsFloat()[i]);
				idx = (uint32_t)((value - minIntensityValue) / (maxIntensityValue - minIntensityValue + 1) * 256);
				break;
			}

			hist[idx]++;
		}

		if (callbackLoadProgressStop != nullptr)
			callbackLoadProgressStop();
	}

	void* GetBuffer()
	{
		itk::simple::PixelIDValueEnum pixelID = GetPixelID();
		switch (pixelID)
		{
		case itk::simple::PixelIDValueEnum::sitkUInt16:
			return(volume->GetBufferAsUInt16());
		case itk::simple::PixelIDValueEnum::sitkInt16:
			return(volume->GetBufferAsInt16());
		case itk::simple::PixelIDValueEnum::sitkUInt8:
			return(volume->GetBufferAsUInt8());
		default:
			throw std::exception("unsupported volume data type encountered...");
		}		
	}

	void PropagateLoadProgress(const std::string caption, float progress) const
	{
		if (callbackOnLoadProgress != nullptr)
			callbackOnLoadProgress(caption, progress);
	}

	std::vector<unsigned int> GetVolumeSize() { return(volume->GetSize()); }

	std::vector<double> GetVoxelSize() { return(volume->GetSpacing()); }

	std::vector<double> GetDirection() { return(volume->GetDirection()); }

	itk::simple::PixelIDValueEnum GetPixelID() { return(volume->GetPixelID()); }
	std::string GetPixelIDTypeAsString() { return(volume->GetPixelIDTypeAsString()); }
};
