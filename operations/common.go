// Copyright (C) 2018, Michael P. Gerlek (Flaxen Consulting)
//
// Portions of this code were derived from the PROJ.4 software
// In keeping with the terms of the PROJ.4 project, this software
// is provided under the MIT-style license in `LICENSE.md` and may
// additionally be subject to the copyrights of the PROJ.4 authors.

package operations

type mode int

const (
	modeNPole mode = 0
	modeSPole      = 1
	modeEquit      = 2
	modeObliq      = 3
)

const tol7 = 1.e-7
const tol10 = 1.0e-10

const eps7 = 1.0e-7
const eps10 = 1.e-10

const FORTPI = 0.78539816339744833
const HALF_PI = 1.570796326794896619
const PI = 3.141592653589793238
const TWO_PI = 6.283185307179586477
const DEL_TOL = 1.0e-14
const MAX_ITER = 20
const D2R = 0.01745329251994329577
