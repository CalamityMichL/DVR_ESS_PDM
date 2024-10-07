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

#include <memory>
#include <fstream>
#include <regex>
#include <format>
#include <span>

#include "SimpleITK.h"

namespace sitk = itk::simple;

class ITF {
public:
	struct SupportPoint {
		unsigned int x;
		unsigned int a, r, g, b; // must be int datatype not uchar8, otherwise stringstream will not correctly convert string later
	};
	
private:	
	std::unique_ptr<sitk::Image> itf_rgba;
	std::unique_ptr<sitk::Image> itf;
	bool applyGammaCorrection;
	float gamma;
	std::vector<SupportPoint> supportPoints;

	std::vector<std::string> split(const std::string str, const std::string regex_str)
	{
		std::regex regexz(regex_str);
		std::vector<std::string> list(std::sregex_token_iterator(str.begin(), str.end(), regexz, -1),
			std::sregex_token_iterator());
		return list;
	}

public:
	ITF()
	{
		itf_rgba = std::make_unique<sitk::Image>(4,256, sitk::sitkUInt8);
		itf = std::make_unique<sitk::Image>(256, 256, sitk::sitkUInt32);
		supportPoints.clear();
		applyGammaCorrection = false;
		gamma = 1.0f;
	}

	void LoadITF(const std::string filepath)
	{
		std::fstream filestream;

		filestream.open(filepath, std::ios::in);

		if (!filestream.is_open())
		{
			throw std::runtime_error("error in ITF::LoadITF: could not load input file " + filepath + "\n");
		}

		supportPoints.clear();

		// first parse header
		std::string line;
		while (getline(filestream, line))
		{
			if (line == "") // an empty line separates the header from the section containing the itf support points
				break;

			std::vector<std::string> tokens = split(line, ":|\n"); // split string with regex, delimiters being ":" and "\n"
			if (tokens[0] == "applyGammaCorrection")
			{
				std::stringstream ss(tokens[1]);
				ss >> applyGammaCorrection;
			}
				
			if (tokens[0] == "gamma")
			{
				std::stringstream ss(tokens[1]);
				ss >> gamma;
			}
		}

		// now parse section containing the itf support points
		while (getline(filestream, line))
		{
			std::vector<std::string> tokens = split(line, ",|\n"); // split string with regex, delimiters being "," and "\n"
			if (tokens.size() == 5) // test if line is a support point definition
			{
				std::stringstream ss;
				SupportPoint supportPoint;
				ss << tokens[0] << ' ' << tokens[1] << ' ' << tokens[2] << ' ' << tokens[3] << ' ' << tokens[4];
				ss >> supportPoint.x >> supportPoint.a >> supportPoint.r >> supportPoint.g >> supportPoint.b;
				supportPoints.push_back(supportPoint);
			}
		}

		filestream.close();

		// and now fill itf_rgba buffer with correct values
		for (size_t i = 0; i < supportPoints.size() - 1; i++)
		{
			SupportPoint p1 = supportPoints[i];
			SupportPoint p2 = supportPoints[i+1];
			SetRangedLinspace(p1.x, p2.x, p1.r, p2.r, p1.g, p2.g, p1.b, p2.b, p1.a, p2.a);
		}
	}

	void SaveITF(const std::string filepath, const std::vector<ITF::SupportPoint> supportPoints, bool applyGamma, float gamma)
	{
		std::string metadata_string = "applyGammaCorrection:" + std::format("{}", static_cast<int>(applyGamma)) + "\n";
		metadata_string += "gamma:" + std::format("{}", gamma) + std::string("\n\n"); // {:.2f}
		
		std::fstream filestream;

		filestream.open(filepath, std::ios::out);
		filestream << metadata_string;
		for (const auto& sp : std::span(supportPoints.begin(), supportPoints.end()))
		{
			filestream <<	std::format("{}", sp.x) + "," +
							std::format("{}", sp.a) + "," +
							std::format("{}", sp.r) + "," +
							std::format("{}", sp.g) + "," +
							std::format("{}", sp.b) + "\n";
		}
		filestream.close();
	}

	void SetSupportPoints(std::vector <SupportPoint> &supportPoints)
	{
		this->supportPoints.clear();
		this->supportPoints = supportPoints;

		// and now fill itf_rgba buffer with correct values
		for (size_t i = 0; i < supportPoints.size() - 1; i++)
		{
			SupportPoint p1 = supportPoints[i];
			SupportPoint p2 = supportPoints[i + 1];
			SetRangedLinspace(p1.x, p2.x, p1.r, p2.r, p1.g, p2.g, p1.b, p2.b, p1.a, p2.a);
		}
	}

	std::vector <SupportPoint> getSupportPoints() const
	{
		std::vector <SupportPoint> retSupportPoints(supportPoints);
		return retSupportPoints;
	}

	void SetRangedLinspace(const unsigned int &lb, const unsigned int &rb, const uint8_t& r_start, const uint8_t& r_end, const uint8_t& g_start, const uint8_t& g_end, const uint8_t& b_start, const uint8_t& b_end, const uint8_t& a_start, const uint8_t& a_end)
	{
		for (unsigned int i = lb; i <= rb; i++)
		{
			float t = (float)(i - lb) / (float)(rb - lb);
			uint8_t rInterp = (uint8_t)std::round(std::lerp((float)r_start, (float)(r_end), t));
			uint8_t gInterp = (uint8_t)std::round(std::lerp((float)g_start, (float)(g_end), t));
			uint8_t bInterp = (uint8_t)std::round(std::lerp((float)b_start, (float)(b_end), t));
			uint8_t aInterp = (uint8_t)std::round(std::lerp((float)a_start, (float)(a_end), t));
			itf_rgba.get()->SetPixelAsUInt8(std::vector<unsigned int> {0, i}, rInterp);
			itf_rgba.get()->SetPixelAsUInt8(std::vector<unsigned int> {1, i}, gInterp);
			itf_rgba.get()->SetPixelAsUInt8(std::vector<unsigned int> {2, i}, bInterp);
			itf_rgba.get()->SetPixelAsUInt8(std::vector<unsigned int> {3, i}, aInterp);
		}
	}

	void SetFromInterpolatedITFs(const ITF &itf1, const ITF &itf2, float t)
	{
		for (unsigned int i = 0; i < 255; i++)
		{
			uint8_t r_interp = (uint8_t)std::round(std::lerp((float)itf1.GetMappedRedValue(i), (float)itf2.GetMappedRedValue(i), t));
			uint8_t g_interp = (uint8_t)std::round(std::lerp((float)itf1.GetMappedGreenValue(i), (float)itf2.GetMappedGreenValue(i), t));
			uint8_t b_interp = (uint8_t)std::round(std::lerp((float)itf1.GetMappedBlueValue(i), (float)itf2.GetMappedBlueValue(i), t));
			uint8_t a_interp = (uint8_t)std::round(std::lerp((float)itf1.GetMappedAlphaValue(i), (float)itf2.GetMappedAlphaValue(i), t));
			itf_rgba->SetPixelAsUInt8(std::vector<unsigned int> {0, i}, r_interp);
			itf_rgba->SetPixelAsUInt8(std::vector<unsigned int> {1, i}, g_interp);
			itf_rgba->SetPixelAsUInt8(std::vector<unsigned int> {2, i}, b_interp);
			itf_rgba->SetPixelAsUInt8(std::vector<unsigned int> {3, i}, a_interp);
		}
		gamma = (float)std::lerp((float)itf1.getGamma(), (float)itf2.getGamma(), t);
		if (gamma > 0.0)
			applyGammaCorrection = true;
		else
			applyGammaCorrection = false;
	}

	uint8_t GetMappedRedValue(uint32_t intensity) const
	{
		return itf_rgba.get()->GetPixelAsUInt8(std::vector<unsigned int> {0, intensity});
	}

	uint8_t GetMappedGreenValue(uint32_t intensity) const
	{
		return itf_rgba.get()->GetPixelAsUInt8(std::vector<unsigned int> {1, intensity});
	}

	uint8_t GetMappedBlueValue(uint32_t intensity) const
	{
		return itf_rgba.get()->GetPixelAsUInt8(std::vector<unsigned int> {2, intensity});
	}

	uint8_t GetMappedAlphaValue(uint32_t intensity) const
	{
		return itf_rgba.get()->GetPixelAsUInt8(std::vector<unsigned int> {3, intensity});
	}

	uint32_t* GetData() const
	{
		for (unsigned int i = 0; i < 255; i++)
		{
			uint32_t value =	((uint32_t)itf_rgba->GetPixelAsUInt8(std::vector<unsigned int> {0, i}) << 24) +
								((uint32_t)itf_rgba->GetPixelAsUInt8(std::vector<unsigned int> {1, i}) << 16) +
								((uint32_t)itf_rgba->GetPixelAsUInt8(std::vector<unsigned int> {2, i}) << 8) +
								((uint32_t)itf_rgba->GetPixelAsUInt8(std::vector<unsigned int> {3, i}));
			itf.get()->SetPixelAsUInt32(std::vector<unsigned int> {i, 0}, value);
			
		}
		return itf.get()->GetBufferAsUInt32();
	}

	bool getApplyGammaCorrection() const {
		return(applyGammaCorrection);
	}

	float getGamma() const {
		return(gamma);
	}
};
