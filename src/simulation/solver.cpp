#include "solver.h"

#include <Eigen/Core>

using Eigen::Vector3f;

// External Force does not changed.

// Function to calculate the derivative of KineticState
KineticState derivative(const KineticState& state)
{
    return KineticState(state.velocity, state.acceleration, Eigen::Vector3f(0, 0, 0));
}

// Function to perform a single Forward Euler step
KineticState forward_euler_step([[maybe_unused]]const KineticState& previous, const KineticState& current)
{
   KineticState future;
    // 计算下一步的速度和位置
  // future.acceleration = current.acceleration;
    future.velocity = current.velocity + (current.acceleration*time_step) ;
    future.position = current.position + (current.velocity * time_step);
    future.acceleration = current.acceleration;
    // 返回下一步的运动学状态
    return future;
    //return current;
} 

// Function to perform a single Runge-Kutta step
KineticState runge_kutta_step([[maybe_unused]] const KineticState& previous,
                              const KineticState& current)
{
    KineticState future;
     // 使用 Runge-Kutta 方法计算下一个时间步长的状态
    Eigen::Vector3f k1_velocity = time_step * current.acceleration;
    Eigen::Vector3f k1_position = time_step * current.velocity;
    Eigen::Vector3f k2_velocity = time_step * (current.acceleration + 0.5f * k1_velocity);
    Eigen::Vector3f k2_position = time_step * (current.velocity + 0.5f * k1_position);
    Eigen::Vector3f k3_velocity = time_step * (current.acceleration + 0.5f * k2_velocity);
    Eigen::Vector3f k3_position = time_step * (current.velocity + 0.5f * k2_position);
    Eigen::Vector3f k4_velocity = time_step * (current.acceleration + k3_velocity);
    Eigen::Vector3f k4_position = time_step * (current.velocity + k3_position);
    
    future.acceleration = current.acceleration;
    future.velocity = current.velocity + (1.0 / 6.0f) * (k1_velocity + 2.0f * k2_velocity + 2.0f * k3_velocity + k4_velocity);
    future.position = current.position + (1.0 / 6.0f) * (k1_position + 2.0f * k2_position + 2.0f * k3_position + k4_position);

    return future;
}

// Function to perform a single Backward Euler step
KineticState backward_euler_step([[maybe_unused]] const KineticState& previous,
                                 const KineticState& current)
{
    KineticState future;
    future.acceleration = current.acceleration;
    future.velocity = current.velocity + (future.acceleration * time_step);
    future.position = current.position + (future.velocity * time_step);
    return future;
}

// Function to perform a single Symplectic Euler step
KineticState symplectic_euler_step(const KineticState& previous, const KineticState& current)
{
    (void)previous;
    KineticState future;
    future.acceleration = current.acceleration;
    future.velocity = current.velocity + (current.acceleration * time_step);
    future.position = current.position + (future.velocity * time_step);
    return future;
}
