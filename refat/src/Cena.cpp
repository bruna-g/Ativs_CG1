#include "../include/Cena.h"

#include "../include/Objeto.h"
#include "../include/Vetor.h"

#include <cmath>

bool Cena::estaNaSombra(const Point& Pi_mod, const Vector& l, float dist_Pi_Luz, Objeto* ignorar) const {
    Vetor origem(Pi_mod.x, Pi_mod.y, Pi_mod.z, 0.0f);
    Vetor dir(l.i, l.j, l.k, 0.0f);

    for (Objeto* obj : objetosSombra) {
        if (obj == nullptr || obj == ignorar) continue;

        if (obj->verificarIntersecao(origem, dir)) {
            float t = static_cast<float>(obj->getDistancia());
            if (t > 1e-4f && std::isfinite(t) && t < dist_Pi_Luz) {
                return true;
            }
        }
    }

    return false;
}
