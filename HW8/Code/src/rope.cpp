#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.
        for(int i = 0; i < num_nodes; ++i)
        {
            float x = float(i) / float(num_nodes);
            Mass *mass = new Mass(start * (1 - x) + end * x, node_mass, false);
            masses.push_back(mass);
        }

        for(int i = 0; i < num_nodes-1; ++i)
        {
            Spring *spring = new Spring(masses[i], masses[i+1], k);
            springs.push_back(spring);
        }

        for (auto &i : pinned_nodes) {
            masses[i]->pinned = true;
        } 
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        const float kd = 0.005;
        for (auto &s : springs)
        {
            Vector2D a(s->m1->position), b(s->m2->position);
            float length = (a - b).norm();
            float k = s->k;
            float rest_length = s->rest_length;
            Vector2D f_a2b = k * (b - a) / length * (length - rest_length);
            s->m1->forces += f_a2b;
            s->m2->forces += -f_a2b;
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position
                m->forces += gravity * m->mass;
                // TODO (Part 2): Add global damping
                m->forces += -kd * m->velocity;
                // Simulate
                Vector2D a = m->forces / m->mass;
                Vector2D v = m->velocity + a * delta_t;
                Vector2D x;
                // Explicit Euler
                //x = m->position + m->velocity * delta_t;
                // Semi-implicit Euler
                x = m->position + v * delta_t;
                m->position = x;
                m->velocity = v;
            }

            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        const float kd = 0.00005;
        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            Vector2D a(s->m1->position), b(s->m2->position);
            float length = (a - b).norm();
            float k = s->k;
            float rest_length = s->rest_length;
            Vector2D f_a2b = k * (b - a) / length * (length - rest_length);
            s->m1->forces += f_a2b;
            s->m2->forces += -f_a2b;
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // TODO (Part 3.1): Set the new position of the rope mass
                // TODO (Part 4): Add global Verlet damping
                m->forces += gravity * m->mass;
                Vector2D a = m->forces / m->mass;

                Vector2D next_position = m->position + (1 - kd) * (m->position - m->last_position) + a * delta_t * delta_t;
                
                m->last_position = m->position;
                m->position = next_position;
            }
            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }
}
