### TO DO

- [ ] Consider spring forces
- [ ] Motility gradient from leader cells to followers
- [ ] Consider no leader identification
- [ ] Initialize with uniform density
- [ ] Number cells in sequence
 
Implemented

- [x] Interaction forces (Spring like adhesion and repulsion; Szabo et al)
- [x] Noise force (Density dependent, mass=1, acceleration is considered, Nirgov)
- [x] Viscek force (Nirgov)
- [x] Boundary detection

To do
- [ ] Proper initialization of positions; semicircular geometry 
- [ ] Incoming cell flux
- [ ] Check density independent noise
- [ ] Check with solely Szabo model 
- [ ] Find boundary cell sequence, and implement extra adhesive interactions 
- [ ] No leader in advance, if leaders emerge when the extra adhesive interactions among the boundary cells are perturbed randomly
- [ ] We may have to consider the model by H Levine (PhysBio 2020; regarding leader cell detection and force laws).

