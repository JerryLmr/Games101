#include <stdexcept>
#include <vector>

#include "rope.h"

namespace CGL {

// Accumulate equal and opposite Hooke forces on each spring's endpoints.
static void applySpringForces(const std::vector<Spring *> &springs)
{
    for (auto *s : springs) {
        const Vector2D displacement = s->m2->position - s->m1->position;
        const double length = displacement.norm();
        if (length == 0.0)
            continue; // A coincident pair has no defined spring direction.
        const Vector2D force = s->k * (length - s->rest_length) * displacement / length;
        s->m1->forces += force;
        s->m2->forces -= force;
    }
}

Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass,
           float k, vector<int> pinned_nodes)
{
    if (num_nodes < 2 || node_mass <= 0 || k < 0)
        throw std::invalid_argument("Rope requires at least two nodes, positive mass and nonnegative stiffness");
    for (int index : pinned_nodes) {
        if (index < 0 || index >= num_nodes)
            throw std::out_of_range("Pinned node index is outside the rope");
    }

    const Vector2D spacing = (end - start) / (num_nodes - 1);
    for (int i = 0; i < num_nodes; ++i) {
        masses.push_back(new Mass(start + i * spacing, node_mass, false));
        if (i > 0)
            springs.push_back(new Spring(masses[i - 1], masses[i], k));
    }
    for (int index : pinned_nodes)
        masses[index]->pinned = true;
}

void Rope::simulateEuler(float delta_t, Vector2D gravity)
{
    applySpringForces(springs);
    const double damping = 0.01;
    for (auto *m : masses) {
        if (!m->pinned) {
            m->forces += m->mass * gravity - damping * m->velocity;
            const Vector2D acceleration = m->forces / m->mass;

            // Explicit Euler: update position using the OLD velocity.
            // m->position += m->velocity * delta_t;
            // m->velocity += acceleration * delta_t;

            // Semi-implicit Euler: update velocity first, then position.
            m->velocity += acceleration * delta_t;
            m->position += m->velocity * delta_t;
        }
        m->forces = Vector2D(0, 0);
    }
}

void Rope::simulateVerlet(float delta_t, Vector2D gravity)
{
    applySpringForces(springs);
    const double damping = 0.00005;
    for (auto *m : masses) {
        if (!m->pinned) {
            const Vector2D previous_position = m->position;
            const Vector2D acceleration = gravity + m->forces / m->mass;
            m->position += (1.0 - damping) * (m->position - m->last_position)
                         + acceleration * delta_t * delta_t;
            m->last_position = previous_position;
        }
        m->forces = Vector2D(0, 0);
    }
}

} // namespace CGL
