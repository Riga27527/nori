#pragma once

#include <nori/common.h>
#include <nori/color.h>
#include <nori/vector.h>
#include <iostream>
#include <opencv2/opencv.hpp>
NORI_NAMESPACE_BEGIN

class Texture{
public:
	Texture() = default;
	
	~Texture() = default;

	Texture(const std::string& filename){
		image_data = cv::imread(filename);
		cv::cvtColor(image_data, image_data, cv::COLOR_BGR2RGB);
		iw = image_data.cols;
		ih = image_data.rows;
		ic = image_data.channels();
	}

	Texture(const Texture& src) = default;

	Texture& operator= (const Texture& src) = default;

	Color3f getPixel(int x, int y) const{
		auto color = image_data.at<cv::Vec3b>(y, x);
		return Color3f(color[0] / 255.0f, color[1] / 255.0f, color[2] / 255.0f);
	}

	Color3f sample2D(float u, float v) const{
		return sampleBilinear(u * iw, (1.0f - v) * ih);
	}

	Color3f sample2D(Point2f uv) const{
		return sample2D(uv[0], uv[1]);
	}

	Color3f sampleBilinear(float x, float y) const;

	int getW() const { return iw; }
	int getH() const { return ih; }
	int getChannel() const { return ic; }

private:
	cv::Mat image_data;
	int iw, ih, ic;
};

NORI_NAMESPACE_END