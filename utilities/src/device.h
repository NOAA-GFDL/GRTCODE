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

#ifndef DEVICE_H_
#define DEVICE_H_

#include "extern.h"


/** @brief Default device ids.*/
enum grtcode_device_ids
{
    HOST_ONLY = -1, /**< Host CPU.*/
    DEFAULT_GPU = 0 /**< Default GPU.*/
};


/** @brief Device.*/
typedef int Device_t;


/** @brief Set the device identifier.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int create_device(Device_t * const device, /**< Device object.*/
                         int const * const id /**< Device identifier.*/
                        );


/** @brief Determine the number of CUDA-enabled GPUs on the system.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int get_num_gpus(int * num_devices, /**< Number of CUDA-enabled devices found.*/
                        int const verbose /**< Verbosity flag.*/
                       );


#endif
