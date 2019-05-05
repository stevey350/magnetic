
#include <iostream>
#include "sysdef.h"
#include "AsyncSerial.h"


using namespace std;
using namespace boost;

#define BUF_SIZE    512

#define MAGIC_NUM1 0xDE
#define MAGIC_NUM2 0xAD
#define MAGIC_NUM3 0xBE
#define MAGIC_NUM4 0xEF

unsigned char buf[BUF_SIZE];
float sensors_data[9][3] = {0};
double valid_sensors_data[9][3] = {0};
short buf_index = 0;
int read_frame_ok = 0;
int geomagnet_removed_ok = 0;

#define LOOP_NUM    20
float geomagnet[SENSOR_NUM][SENSOR_DIM] = {0};
int cnt_num = 0;

void readCallback(const char *data_packet, size_t size)
{
    if( (buf_index + size) > BUF_SIZE)
    {
        cout << "buf overflowed." <<endl;
        buf_index = 0;
    }
    memcpy(buf+buf_index, data_packet, size);
    buf_index += size;

    // is end ?
//    cout << "data: ";
//    for(int i=0; i<size; i++)
//        cout << hex << int(data_packet[i]);
//    cout << endl;
    if( (buf[buf_index-1] == MAGIC_NUM4) && (buf[buf_index-2] == MAGIC_NUM3) \
            && (buf[buf_index-3] == MAGIC_NUM2) && (buf[buf_index-4] == MAGIC_NUM1) )
    {
        if(buf_index == (sizeof(sensors_data)+4))
        {
            memcpy(sensors_data, buf, buf_index-4);
//            for(int j=0; j<9; j++)
//                cout << j << "- x: " << sensors_data[j][0] << " y: " << sensors_data[j][1] << " z: " << sensors_data[j][2] << endl;

            // remove geomagnet
            if(cnt_num < LOOP_NUM)
            {
                cnt_num++;
                for(int i=0; i<SENSOR_NUM; i++)
                {
                    geomagnet[i][0] += sensors_data[i][0];
                    geomagnet[i][1] += sensors_data[i][1];
                    geomagnet[i][2] += sensors_data[i][2];
                }

                if(LOOP_NUM == cnt_num)
                {
                    for(int i=0; i<SENSOR_NUM; i++)
                    {
                        geomagnet[i][0] /= LOOP_NUM;
                        geomagnet[i][1] /= LOOP_NUM;
                        geomagnet[i][2] /= LOOP_NUM;
                    }
                    cout << "geomagnet has been removed" << endl;
                    geomagnet_removed_ok = 1;
                }
                return;
            }
            for(int i=0; i<SENSOR_NUM; i++)
            {
                sensors_data[i][0] = sensors_data[i][0] - geomagnet[i][0];
                sensors_data[i][1] = sensors_data[i][1] - geomagnet[i][1];
                sensors_data[i][2] = sensors_data[i][2] - geomagnet[i][2];
            }

            // correct the sensors direction
            for (int i=0; i<SENSOR_NUM; i++)
            {
                valid_sensors_data[i][0] = sensors_data[i][0]*sensors_dir[i*SENSOR_DIM+0][0]+
                                     sensors_data[i][1]*sensors_dir[i*SENSOR_DIM+0][1]+
                                     sensors_data[i][2]*sensors_dir[i*SENSOR_DIM+0][2];

                valid_sensors_data[i][1] = sensors_data[i][0]*sensors_dir[i*SENSOR_DIM+1][0]+
                                     sensors_data[i][1]*sensors_dir[i*SENSOR_DIM+1][1]+
                                     sensors_data[i][2]*sensors_dir[i*SENSOR_DIM+1][2];

                valid_sensors_data[i][2] = sensors_data[i][0]*sensors_dir[i*SENSOR_DIM+2][0]+
                                     sensors_data[i][1]*sensors_dir[i*SENSOR_DIM+2][1]+
                                     sensors_data[i][2]*sensors_dir[i*SENSOR_DIM+2][2];

                valid_sensors_data[i][0] = valid_sensors_data[i][0] * 1e-6;     // The unit of uT converts to T
                valid_sensors_data[i][1] = valid_sensors_data[i][1] * 1e-6;
                valid_sensors_data[i][2] = valid_sensors_data[i][2] * 1e-6;
            }

        //    cout << "is ok" << endl;
            read_frame_ok = 1;
        }

        buf_index = 0;
    }
}

