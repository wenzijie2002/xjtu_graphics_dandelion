#include "bvh.h"

#include <cassert>

#include <optional>

#include <Eigen/Geometry>
#include "formatter.hpp"
#include <spdlog/spdlog.h>

#include "math.hpp"

using Eigen::Vector3f;
using std::optional;
using std::vector;

BVHNode::BVHNode() : left(nullptr), right(nullptr), face_idx(0)
{
}

BVH::BVH(const GL::Mesh& mesh) : root(nullptr), mesh(mesh)
{
}

// 建立bvh，将需要建立BVH的图元索引初始化
void BVH::build()
{
	if (mesh.faces.count() == 0) {
		root = nullptr;
		return;
	}

	primitives.resize(mesh.faces.count());
	for (size_t i = 0; i < mesh.faces.count(); i++) primitives[i] = i;

	root = recursively_build(primitives);
	return;
}
// 删除bvh
void BVH::recursively_delete(BVHNode* node)
{
	if (node == nullptr)
		return;
	recursively_delete(node->left);
	recursively_delete(node->right);
	delete node;
	node = nullptr;
}
// 统计BVH树建立的节点个数
size_t BVH::count_nodes(BVHNode* node)
{
	if (node == nullptr)
		return 0;
	else
		return count_nodes(node->left) + count_nodes(node->right) + 1;
}
// 递归建立BVH
BVHNode* BVH::recursively_build(vector<size_t> faces_idx)
{
	BVHNode* node = new BVHNode();

	AABB aabb;
	for (size_t i = 0; i < faces_idx.size(); i++) {
		aabb = union_AABB(aabb, get_aabb(mesh, faces_idx[i]));

	}
	if (faces_idx.size() == 1) {
		// 如果只有一个面，返回叶子节点
		node->aabb = aabb;
		node->face_idx = faces_idx[0];
		node->left = nullptr;
		node->right = nullptr;
		return node;
	}
	else if (faces_idx.size() == 2) {
	 // 如果有两个面，分别递归构建左右节点
		node->left = recursively_build(vector<size_t>{faces_idx[0]});
		node->right = recursively_build(vector<size_t>{faces_idx[1]});
		// 更新当前节点的AABB
		node->aabb = union_AABB(node->left->aabb, node->right->aabb);
		return node;
	}
	else {
	 // 选择x、y、z中最长的一维
		int dim = aabb.max_extent();
		// 沿最长维度对面进行排序
		sort(faces_idx.begin(), faces_idx.end(), [dim, this](size_t a, size_t b) {
			return get_aabb(mesh, a).centroid()[dim] < get_aabb(mesh, b).centroid()[dim];
			});

		// 沿最长维度将图元分成两部分
		std::vector<size_t> left_faces_idx(faces_idx.begin(), faces_idx.begin() + faces_idx.size() / 2);
		std::vector<size_t> right_faces_idx(faces_idx.begin() + faces_idx.size() / 2, faces_idx.end());
		// 递归构建左右节点
		node->left = recursively_build(left_faces_idx);
		node->right = recursively_build(right_faces_idx);
		// 更新当前节点的AABB
		node->aabb = union_AABB(node->left->aabb, node->right->aabb);
		return node;
	}
}
optional<Intersection> BVH::intersect(const Ray& ray, [[maybe_unused]] const GL::Mesh& mesh,
	const Eigen::Matrix4f obj_model)
{
	model = obj_model;
	optional<Intersection> isect;
	if (!root) {
		isect = std::nullopt;
		return isect;
	}
	isect = ray_node_intersect(root, ray);
	// 将求交结果从模型坐标系变换到世界坐标系下
	if (isect) {
		Vector3f world_intersection =
			(model * (isect->barycentric_coord).homogeneous()).head<3>();
		isect->t = (world_intersection - ray.origin).norm();
		isect->normal =
			(model * Eigen::Vector4f(isect->normal.x(), isect->normal.y(), isect->normal.z(), 0.0f))
			.head<3>().normalized();
	}
	return isect;
}
optional<Intersection> BVH::ray_node_intersect(BVHNode* node, const Ray& ray) const
{
	
	optional<Intersection> isect=std::nullopt;
	
	// 将射线变换到模型坐标系下
	Ray trans_ray; 
	Eigen::Matrix4f model_inv = model.inverse();
	trans_ray.origin =
		(model_inv * ray.origin.homogeneous()).hnormalized();
	trans_ray.direction = (model_inv * Eigen::Vector4f(ray.direction.x(), ray.direction.y(),
		ray.direction.z(), 0.0f)).head(3);

	Eigen::Vector3f inv_dir = trans_ray.direction.cwiseInverse();
	std::array<int, 3> dir_is_neg;
	dir_is_neg[0] = (trans_ray.direction.x() < 0) ? 0 : 1;
	dir_is_neg[1] = (trans_ray.direction.y() < 0) ? 0 : 1;
	dir_is_neg[2] = (trans_ray.direction.z() < 0) ? 0 : 1;

	// 检查射线与当前节点的包围盒是否相交
	if (node->aabb.intersect(trans_ray, inv_dir, dir_is_neg)) {
		if (node->left == nullptr && node->right == nullptr) {
			// 如果是叶节点，直接计算与面片的交点
			isect = ray_triangle_intersect(trans_ray, mesh, node->face_idx);
		}
		else {
			// 非叶子节点，递归地检查左右子节点
			optional<Intersection> left_isect, right_isect;
			left_isect = ray_node_intersect(node->left, ray);
			right_isect = ray_node_intersect(node->right, ray);
			isect= left_isect.has_value() ? left_isect : right_isect;
		}
	}
	return isect;
}


