#include <chemfiles.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef DEBUG
#define DPRINT(...)                                                            \
  do {                                                                         \
    fprintf(stderr, __VA_ARGS__);                                              \
  } while (0)
#else
#define DPRINT(...)
#endif

void readTrajectory(CHFL_TRAJECTORY *file, unsigned long nblocksteps,
                    size_t natoms, float *block_pos, float *block_vel,
                    float *block_box) {

  // for reading of frame
  CHFL_FRAME *frame;
  CHFL_CELL *cell;
  chfl_vector3d *r = NULL;
  chfl_vector3d *v = NULL;
  uint64_t natoms_traj = 0;
  chfl_vector3d box = {0, 0, 0};

  DPRINT("start reading frame\n");
  for (unsigned long t = 0; t < nblocksteps; t++) {
    frame = chfl_frame();
    chfl_trajectory_read(file, frame);
    chfl_frame_positions(frame, &r, &natoms_traj);
    chfl_frame_velocities(frame, &v, &natoms_traj);
    cell = chfl_cell_from_frame(frame);
    chfl_cell_lengths(cell, box);

    DPRINT("There are %lu (ignoring %lu) atoms at step %lu. My box is: %f %f "
           "%f \n",
           natoms_traj, natoms_traj - natoms, t, box[0], box[1], box[2]);

    block_box[3 * t + 0] = box[0] / 10.0;
    block_box[3 * t + 1] = box[1] / 10.0;
    block_box[3 * t + 2] = box[2] / 10.0;
    for (size_t j = 0; j < natoms; j++) {
      block_pos[3 * natoms * t + 3 * j + 0] = r[j][0] / 10.0;
      block_pos[3 * natoms * t + 3 * j + 1] = r[j][1] / 10.0;
      block_pos[3 * natoms * t + 3 * j + 2] = r[j][2] / 10.0;
      block_vel[3 * natoms * t + 3 * j + 0] = v[j][0] / 10.0;
      block_vel[3 * natoms * t + 3 * j + 1] = v[j][1] / 10.0;
      block_vel[3 * natoms * t + 3 * j + 2] = v[j][2] / 10.0;
    }
    // free stuff
    chfl_free(cell);
    chfl_free(frame);
  }
}
