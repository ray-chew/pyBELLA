Terrain phase 2: pressure gradients map to physical space via the terrain
gradient matrix A (slope correction of the horizontal rows, 1/J on the
vertical) in both the explicit forward step and the implicit pressure
correction; the explicit π update divides the J-weighted divergence by the
node Jacobian. Bypassed entirely without terrain.
