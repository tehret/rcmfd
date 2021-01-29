/*
 * Copyright (c) 2019, Thibaud Ehret <ehret.thibaud@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU Affero General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef WRAPPER_RCMFD_H
#define WRAPPER_RCMFD_H

bool perform_matching_py(int channels, float* image, int w, int h, int ps, float tau, bool automatic, float* output, bool verbose);

#endif // WRAPPER_RCMFD_H
