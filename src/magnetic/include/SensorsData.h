
#ifndef __SENSORSDATA_H_
#define __SENSORSDATA_H_

extern double valid_sensors_data[9][3];
extern int read_frame_ok;
extern int geomagnet_removed_ok;

void readCallback(const char *data_packet, size_t size);



#endif
