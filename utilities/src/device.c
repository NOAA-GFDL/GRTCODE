/* GRTCODE is a GPU-able Radiative Transfer Code
 * Copyright (C) 2016  Garrett Wright
 * Modified in 2019 by Raymond Menzel
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; version 2.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "debug.h"
#include "device.h"
#include "extern.h"
#include "return_codes.h"


/*Determine the number of CUDA-enabled GPUs on the system.*/
EXTERN int get_num_gpus(int * num_devices, int const verbose)
{
    not_null(num_devices);
#ifdef __NVCC__
    gpu_catch(cudaGetDeviceCount(num_devices));
    if (verbose)
    {
        char const *mesg = "Found %d GPU devices:";
        log_mesg(mesg, *num_devices);
        int i;
        for (i=0; i<(*num_devices); ++i)
        {
            cudaDeviceProp prop;
            gpu_catch(cudaGetDeviceProperties(&prop, i));
            mesg = "\tDevice #%d: %s";
            log_mesg(mesg, i, prop.name);
        }
    }
#else
    (void)(verbose); /*To quiet compiler warning.*/
    *num_devices = 0;
#endif
    return GRTCODE_SUCCESS;
}


/*Set the device identifier.*/
EXTERN int create_device(Device_t * const device, int const * const id)
{
    not_null(device);
    int num_devices;
    catch(get_num_gpus(&num_devices, 1));
    if (id != NULL)
    {
        if (*id != HOST_ONLY)
        {
            in_range(*id, 0, num_devices);
        }
        *device = *id;
    }
    else if (num_devices > 0)
    {
        *device = DEFAULT_GPU;
    }
    else
    {
        *device = HOST_ONLY;
    }
    return GRTCODE_SUCCESS;
}
