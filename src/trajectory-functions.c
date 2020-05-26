#include <chemfiles.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void get_frame_pos_box(CHFL_FRAME *frame, size_t natoms, float *pos, float *box) {

    // for reading of frame
    CHFL_CELL *cell;
    chfl_vector3d *r = NULL;
    uint64_t natoms_traj = 0;
    chfl_vector3d box_temp = {0, 0, 0};

    chfl_frame_positions(frame, &r, &natoms_traj);
    cell = chfl_cell_from_frame(frame);
    chfl_cell_lengths(cell, box_temp);

    box[0] = box_temp[0] / 10.0;
    box[1] = box_temp[1] / 10.0;
    box[2] = box_temp[2] / 10.0;
    for (size_t j = 0; j < natoms; j++) {
        pos[3 * j + 0] = r[j][0] / 10.0;
        pos[3 * j + 1] = r[j][1] / 10.0;
        pos[3 * j + 2] = r[j][2] / 10.0;
    }
    // free stuff
    chfl_free(cell);
}

void get_traj_pos_vel_box(CHFL_TRAJECTORY *file, unsigned long nblocksteps,
                         size_t natoms, float *block_pos, float *block_vel,
                         float *block_box) {

  // for reading of frame
  CHFL_FRAME *frame;
  CHFL_CELL *cell;
  chfl_vector3d *r = NULL;
  chfl_vector3d *v = NULL;
  uint64_t natoms_traj = 0;
  chfl_vector3d box = {0, 0, 0};

  for (unsigned long t = 0; t < nblocksteps; t++) {
    frame = chfl_frame();
    chfl_trajectory_read(file, frame);
    chfl_frame_positions(frame, &r, &natoms_traj);
    chfl_frame_velocities(frame, &v, &natoms_traj);
    cell = chfl_cell_from_frame(frame);
    chfl_cell_lengths(cell, box);

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

void check_frame_natoms(
        CHFL_FRAME *frame,
        size_t natoms)
{
    uint64_t natoms_traj = 0;
    chfl_vector3d *positions = NULL;
    chfl_frame_positions(frame, &positions, &natoms_traj);
    if (natoms_traj < natoms) {
        fprintf(stderr, "ERROR: The topology you give has more atoms than first "
                "frame of the trajectory/refconf\n");
        exit(1);
    } else if (natoms_traj > natoms) {
        fprintf(stderr, "WARNING: The topology you give has less atoms than first "
                "frame of the trajectory/refconf\n");
        fprintf(stderr, "         Some atoms are ignored in every frame\n");
    }
}

void check_frame_velocities(
        CHFL_FRAME *frame)
{
    bool has_velocities = false;
    chfl_frame_has_velocities(frame, &has_velocities);
    if (!has_velocities) {
        fprintf(stderr, "ERROR: No velocities in trajectory.\n");
        exit(1);
    }
}

void check_frame_orthorombic_box(
        CHFL_FRAME *frame,
        bool no_pbc)
{
    CHFL_CELL *cell = chfl_cell_from_frame(frame);
    chfl_cellshape shape;
    chfl_cell_shape(cell, &shape);
    if ((no_pbc == false) && (shape != CHFL_CELL_ORTHORHOMBIC)) {
        fprintf(stderr,
                "ERROR: can not do recombination on non orthorombic box.\n");
        exit(1);
    }
    chfl_free(cell);
}

float get_traj_framelength(
        CHFL_TRAJECTORY *file,
        CHFL_FRAME *frame)
{
    float framelength;
    double time0, time1;
    int result0, result1;
    // get time0
    CHFL_PROPERTY *property = chfl_frame_get_property(frame, "time");
    result0 = chfl_property_get_double(property, &time0);
    chfl_free(property);
    // get time1
    chfl_trajectory_read(file, frame);
    property = chfl_frame_get_property(frame, "time");
    result1 = chfl_property_get_double(property, &time1);
    chfl_free(property);
    framelength = (float)(time1 - time0);
    if ((framelength == 0.0) || (result0 != CHFL_SUCCESS) ||
        (result1 != CHFL_SUCCESS)) {
      fprintf(stderr,
              "ERROR: Reading framelength from trajectory failed. You can "
              "provide the framelength with command line arguments.\n");
      exit(1);
    }
    return framelength;
}
