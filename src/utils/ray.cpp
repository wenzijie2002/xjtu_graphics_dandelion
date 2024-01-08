#include "ray.h"

#include <cmath>
#include <array>

#include <Eigen/Dense>
#include <spdlog/spdlog.h>

#include "../utils/math.hpp"

using Eigen::Matrix3f;
using Eigen::Matrix4f;
using Eigen::Vector2f;
using Eigen::Vector3f;
using Eigen::Vector4f;
using std::numeric_limits;
using std::optional;
using std::size_t;

constexpr float infinity = 1e5f;
constexpr float eps = 1e-5f;

Intersection::Intersection() : t(numeric_limits<float>::infinity()), face_index(0)
{
}

Ray generate_ray(int width, int height, int x, int y, Camera& camera, float depth)
{
	// these lines below are just for compiling and can be deleted
	(void)width;
	(void)height;
	(void)x;
	(void)y;
	(void)depth;
	// these lines above are just for compiling and can be deleted


	// The ratio between the specified plane (width x height)'s depth and the image plane's depth.

	// Transfer the view-space position to world space.
	Vector3f world_pos;
	return { camera.position, (world_pos - camera.position).normalized() };
}

optional<Intersection> ray_triangle_intersect(const Ray& ray, const GL::Mesh& mesh, size_t index)
{

	Intersection result;
	float min_t = infinity;

	// 获取三角形的顶点
	const auto& face = mesh.face(index);
	Vector3f a = mesh.vertex(face[0]);
	Vector3f b = mesh.vertex(face[1]);
	Vector3f c = mesh.vertex(face[2]);

	Vector3f ab = b - a; // 计算ab向量
	Vector3f ac = c - a; // 计算ac向量
	Vector3f normal = ab.cross(ac).normalized(); // 计算面片法向量并归一化

	Matrix3f A;
	A.col(0) = ray.direction;
	A.col(1) = a - b;
	A.col(2) = a - c;

	//判断是否平行
	if (abs(A.determinant()) <= eps)
	{
		return std::nullopt;
	}

	Vector3f x = A.inverse()*(a - ray.origin); // 解方程Ax = (a - origin)，得到交点参数
	float beta = x[1];
	float gamma = x[2];
	float t = x[0];
	float alpha = 1.0f - beta - gamma; // 计算重心坐标alpha

	// 判断交点是否在三角形内部
	if (beta >= 0 && gamma >= 0 && alpha >= 0 && t > 0 && t < min_t) {
		min_t = t; // 更新最小的t值
		result.t = t; // 更新交点的t值
		result.barycentric_coord = ray.origin + t * ray.direction;// 计算交点坐标
		result.normal = normal; // 更新交点的法向量
	}
	if (min_t < infinity - eps) {
		return result; // 返回最终计算得到的交点
	}
	else {
		return std::nullopt; // 如果没有找到交点，则返回空值
	}
}

optional<Intersection> naive_intersect(const Ray& ray, const GL::Mesh& mesh, const Matrix4f model)
{
	Intersection result;
	float min_t = infinity;

	for (size_t i = 0; i < mesh.faces.count(); ++i) {
		const auto& face = mesh.face(i);

		// 获取三角形的顶点
		Vector3f v0 = mesh.vertex(face[0]);
		Vector3f v1 = mesh.vertex(face[1]);
		Vector3f v2 = mesh.vertex(face[2]);


		const Vector3f& a = (model * v0.homogeneous()).hnormalized();
		const Vector3f& b = (model * v1.homogeneous()).hnormalized();
		const Vector3f& c = (model * v2.homogeneous()).hnormalized();

		Vector3f ab = b - a; // 计算ab向量
		Vector3f ac = c - a; // 计算ac向量
		Vector3f normal = ab.cross(ac).normalized(); // 计算面片法向量并归一化

		Matrix3f A;
		A.col(0) = ray.direction;
		A.col(1) = a - b;
		A.col(2) = a - c;

		Vector3f x = A.partialPivLu().solve(a - ray.origin); // 解方程Ax = (a - origin)，得到交点参数
		float beta = x[1];
		float gamma = x[2];
		float t = x[0];
		//std::cout << t << "\n";
		float alpha = 1.0f - beta - gamma; // 计算重心坐标alpha

		// 判断交点是否在三角形内部
		if (A.determinant() != 0 && beta >= 0 && gamma >= 0 && alpha >= 0 && t > 0 && t < min_t) {
			min_t = t; // 更新最小的t值
			result.t = t; // 更新交点的t值
			result.barycentric_coord = ray.origin + t * ray.direction; // 计算交点的世界坐标
			result.normal = normal; // 更新交点的法向量
		}
	}

	if (min_t < infinity - eps) {
		return result; // 返回最终计算得到的交点
	}
	else {
		return std::nullopt; // 如果没有找到交点，则返回空值
	}
}

