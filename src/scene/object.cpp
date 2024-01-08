#include "object.h"

#include <array>
#include <optional>

#ifdef _WIN32
#include <Windows.h>
#endif
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <fmt/format.h>

#include "../utils/math.hpp"
#include "../utils/ray.h"
#include "../simulation/solver.h"
#include "../utils/logger.h"

using Eigen::Matrix4f;
using Eigen::Quaternionf;
using Eigen::Vector3f;
using std::array;
using std::make_unique;
using std::optional;
using std::string;
using std::vector;

bool Object::BVH_for_collision = false;
size_t Object::next_available_id = 0;
std::function<KineticState(const KineticState&, const KineticState&)> Object::step =
forward_euler_step;

Object::Object(const string& object_name)
	: name(object_name), center(0.0f, 0.0f, 0.0f), scaling(1.0f, 1.0f, 1.0f),
	rotation(1.0f, 0.0f, 0.0f, 0.0f), velocity(0.0f, 0.0f, 0.0f), force(0.0f, 0.0f, 0.0f),
	mass(1.0f), BVH_boxes("BVH", GL::Mesh::highlight_wireframe_color)
{
	visible = true;
	modified = false;
	id = next_available_id;
	++next_available_id;
	bvh = make_unique<BVH>(mesh);
	const string logger_name = fmt::format("{} (Object ID: {})", name, id);
	logger = get_logger(logger_name);
}

Matrix4f Object::model()
{
	 // 创建平移矩阵
	Eigen::Matrix4f translate_mat = Eigen::Matrix4f::Identity();
	translate_mat.block<3, 1>(0, 3) = center;

	// 创建旋转矩阵
	Eigen::Matrix4f rotate_mat = Eigen::Matrix4f::Identity();
	Eigen::Matrix4f rotation_x_matrix = Eigen::Matrix4f::Identity();
	Eigen::Matrix4f rotation_y_matrix = Eigen::Matrix4f::Identity();
	Eigen::Matrix4f rotation_z_matrix = Eigen::Matrix4f::Identity();
	const Quaternionf& r = rotation;
	auto [x_angle, y_angle, z_angle] = quaternion_to_ZYX_euler(r.w(), r.x(), r.y(), r.z());
	double x = x_angle * pi<double>() / 180.0; double y = y_angle * pi<double>() / 180.0; double z = z_angle * pi<double>() / 180.0;
	rotation_x_matrix(1, 1) = cos(x); rotation_x_matrix(2, 2) = cos(x); //x
	rotation_x_matrix(1, 2) = -sin(x); rotation_x_matrix(2, 1) = sin(x);
	rotation_y_matrix(0, 0) = cos(y); rotation_y_matrix(2, 2) = cos(y);//y
	rotation_y_matrix(0, 2) = sin(y); rotation_y_matrix(2, 0) = -sin(y);
	rotation_z_matrix(0, 0) = cos(z); rotation_z_matrix(1, 1) = cos(z);//z
	rotation_z_matrix(0, 1) = -sin(z); rotation_z_matrix(1, 0) = sin(z);
	rotate_mat = rotation_x_matrix * rotation_y_matrix * rotation_z_matrix;
	//rotate_mat.block<3, 3>(0, 0) = rotation.toRotationMatrix();

	// 创建缩放矩阵
	Eigen::Matrix4f scale_mat = Eigen::Matrix4f::Identity();
	scale_mat(0, 0) = scaling.x();
	scale_mat(1, 1) = scaling.y();
	scale_mat(2, 2) = scaling.z();

	// 计算最终的变换矩阵
	Eigen::Matrix4f transformation_matrix = translate_mat * rotate_mat * scale_mat;

	return transformation_matrix;
}

void Object::update(vector<Object*>& all_objects)
{
	// 首先调用 step 函数计下一步该物体的运动学状态。
	KineticState current_state{ center, velocity, force / mass };
	KineticState next_state = step(prev_state, current_state);
	// 将物体的位置移动到下一步状态处，但暂时不要修改物体的速度。
	center = next_state.position;
	// 遍历 all_objects，检查该物体在下一步状态的位置处是否会与其他物体发生碰撞。
	if (BVH_for_collision == false) {
		for (auto object : all_objects) {
			if (object != this) {
				// 检测该物体与另一物体是否碰撞的方法是：
				// 遍历该物体的每一条边，构造与边重合的射线去和另一物体求交，如果求交结果非空、
				// 相交处也在这条边的两个端点之间，那么该物体与另一物体发生碰撞。
				// 请时刻注意：物体 mesh 顶点的坐标都在模型坐标系下，你需要先将其变换到世界坐标系。
				for (size_t i = 0; i < mesh.edges.count(); ++i) {
					array<size_t, 2> v_indices = mesh.edge(i);
					// v_indices 中是这条边两个端点的索引，以这两个索引为参数调用 GL::Mesh::vertex
					// 方法可以获得它们的坐标，进而用于构造射线。
					Vector3f p1 = (model() * mesh.vertex(v_indices[0]).homogeneous()).hnormalized();
					Vector3f p2 = (model() * mesh.vertex(v_indices[1]).homogeneous()).hnormalized();
					Ray edge_ray{ p1, p2 - p1 };
					// 检测射线与另一物体是否相交
					auto intersection = naive_intersect(edge_ray, object->mesh, object->model());
					//else{auto intersection=bvh->ray_node_intersect(object->bvh->root,)}
					if (intersection && intersection->t > 0 && intersection->t < 1) {
						// 如果相交，则计算碰撞后的速度
						float other_mass = object->mass;
						Vector3f other_velocity = object->velocity;

						Vector3f relative_velocity = velocity - other_velocity;
						auto normal = intersection->normal;
						auto impulse_magnitude = 2.0f / (1.0f / mass + 1.0f / other_mass) * relative_velocity.dot(normal);

						Vector3f impulse = impulse_magnitude * normal;

						velocity -= impulse / mass;
						object->velocity += impulse / other_mass;

						// 将下一步状态设为当前状态，以避免重复碰撞
						center = current_state.position;
						return;
					}
				}
			}
		}
	}
	else {//使用bvh
		// 遍历 all_objects，检查该物体在下一步状态的位置处是否会与其他物体发生碰撞。
		for (auto object : all_objects) {
			if (object != this) {
				// 检测该物体与另一物体是否碰撞的方法是使用BVH树进行加速碰撞检测
				for (size_t i = 0; i < mesh.edges.count(); ++i) {
					array<size_t, 2> v_indices = mesh.edge(i);
					// v_indices 中是这条边两个端点的索引，以这两个索引为参数调用 GL::Mesh::vertex
					// 方法可以获得它们的坐标，进而用于构造射线。
					Vector3f p1 = (model() * mesh.vertex(v_indices[0]).homogeneous()).hnormalized();
					Vector3f p2 = (model() * mesh.vertex(v_indices[1]).homogeneous()).hnormalized();
					Ray edge_ray{ p1, (p2 - p1)};
					auto intersection = object->bvh->intersect(edge_ray, object->bvh->mesh, object->model());
					if (intersection && intersection->t > 0 && intersection->t < 1) {
						// 如果发生碰撞，处理碰撞逻辑
						// 如果相交，则计算碰撞后的速度
						float other_mass = object->mass;
						Vector3f other_velocity = object->velocity;

						Vector3f relative_velocity = velocity - other_velocity;
						auto normal = intersection->normal;
						auto impulse_magnitude = 2.0f / (1.0f / mass + 1.0f / other_mass) * relative_velocity.dot(normal);

						Vector3f impulse = impulse_magnitude * normal;

						velocity -= impulse / mass;
						object->velocity += impulse / other_mass;
						// 将下一步状态设为当前状态，以避免重复碰撞
						center = current_state.position;
						return;
					}
				}
			}
		}
	}
	// 将上一步状态赋值为当前状态，并将物体更新到下一步状态。
	prev_state = current_state;
	velocity = next_state.velocity;
	force = next_state.acceleration * mass;
}

void Object::render(const Shader& shader, WorkingMode mode, bool selected)
{
	if (modified) {
		mesh.VAO.bind();
		mesh.vertices.to_gpu();
		mesh.normals.to_gpu();
		mesh.edges.to_gpu();
		mesh.edges.release();
		mesh.faces.to_gpu();
		mesh.faces.release();
		mesh.VAO.release();
	}
	modified = false;
	// Render faces anyway.
	unsigned int element_flags = GL::Mesh::faces_flag;
	if (mode == WorkingMode::MODEL) {
		// For *Model* mode, only the selected object is rendered at the center in the world.
		// So the model transform is the identity matrix.
		shader.set_uniform("model", I4f);
		shader.set_uniform("normal_transform", I4f);
		element_flags |= GL::Mesh::vertices_flag;
		element_flags |= GL::Mesh::edges_flag;
	}
	else {
		Matrix4f model = this->model();
		shader.set_uniform("model", model);
		shader.set_uniform("normal_transform", (Matrix4f)(model.inverse().transpose()));
	}
	// Render edges of the selected object for modes with picking enabled.
	if (check_picking_enabled(mode) && selected) {
		element_flags |= GL::Mesh::edges_flag;
	}
	mesh.render(shader, element_flags);
}

void Object::rebuild_BVH()
{
	bvh->recursively_delete(bvh->root);
	bvh->build();
	BVH_boxes.clear();
	refresh_BVH_boxes(bvh->root);
	BVH_boxes.to_gpu();
}

void Object::refresh_BVH_boxes(BVHNode* node)
{
	if (node == nullptr) {
		return;
	}
	BVH_boxes.add_AABB(node->aabb.p_min, node->aabb.p_max);
	refresh_BVH_boxes(node->left);
	refresh_BVH_boxes(node->right);
}
