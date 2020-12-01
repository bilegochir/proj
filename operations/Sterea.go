package operations

import (
	"math"

	"github.com/bilegochir/proj/core"
)

func init() {
	core.RegisterConvertLPToXY("sterea",
		"Oblique Stereographic Alternative",
		"asd",
		NewSterea,
	)
}

// Sterea implements core.IOperation and core.ConvertLPToXY
// Lam0, Phi0     float64 /* central meridian, parallel */
type Sterea struct {
	core.Operation
	// Gauss
	lat0 float64
	sphi float64
	cphi float64
	rc float64
	phic0 float64
	C float64
	ratexp float64
	K float64

	// Sterea
	sinc0 float64
	cosc0 float64
	R2 float64
    long0 float64
}

// NewSterea creates a new Plate Carree system
func NewSterea(system *core.System, desc *core.OperationDescription) (core.IConvertLPToXY, error) {
	op := &Sterea{}
	op.System = system

	err := op.stereaSetup(system)
	if err != nil {
		return nil, err
	}
	return op, nil
}

// Forward goes forewards
func (op *Sterea) Forward(lp *core.CoordLP) (*core.CoordXY, error) {
	xy := &core.CoordXY{X: 0.0, Y: 0.0}
	P := op.System
	PE := op.System.Ellipsoid

	xyz := &core.CoordLPZ{Lam: lp.Lam, Phi: lp.Phi, Z: 0.0}
	xyz = op.geodetic_to_geocentric(xyz)
	xyz = op.geocentric_from_wgs84(xyz, op.System.DatumParams)
	xyz = op.geocentric_to_geodetic(xyz, op.System)

	xyz.Lam = op.adjust_lon(xyz.Lam - op.long0)

	lon := xyz.Lam
	lat := xyz.Phi
	xy.Y = 2.0 * math.Atan(op.K * math.Pow(math.Tan(0.5 * lat + FORTPI), op.C) * op.srat(PE.E * math.Sin(lat), op.ratexp)) - HALF_PI
	xy.X = op.C * lon

	sinc := math.Sin(xy.Y)
	cosc := math.Cos(xy.Y)
	cosl := math.Cos(xy.X)
	k := P.K0 * op.R2 / (1.0 + op.sinc0 * sinc + op.cosc0 * cosc * cosl)

	xy.X = k * cosc * math.Sin(xy.X)
	xy.Y = k * (op.cosc0 * sinc - op.sinc0 * cosc * cosl)

	xy.X = 6377397.155 * xy.X + P.X0
	xy.Y = 6377397.155 * xy.Y + P.Y0

	return xy, nil
}

// Inverse goes backwards
func (op *Sterea) Inverse(xy *core.CoordXY) (*core.CoordLP, error) {
	lp := &core.CoordLP{Lam: 0.0, Phi: 0.0}
	P := op.System
	PE := op.System.Ellipsoid

	lon := lp.Lam / op.C
	lat := lp.Phi
	num := math.Pow(math.Tan(0.5 * lat + FORTPI) / op.K, 1.0 / op.C)

	for i := MAX_ITER; i > 0; i -= 20 {
		lat := 2.0 * math.Atan(num * op.srat(PE.E * math.Sin(lp.Phi), -0.5 * PE.E)) - HALF_PI

		if math.Abs(lp.Lam - lp.Phi) < DEL_TOL {
			break
		}

		lp.Lam = lat
	}

	lp.Lam = lon
	lp.Phi = lat

	lon = (lp.Lam - P.X0) / 6377397.155
	lp.Phi = (lp.Phi - P.Y0) / 6377397.155

	lp.Lam = lp.Lam / P.K0
	lp.Phi = lp.Phi / P.K0

	rho := math.Sqrt(lp.Lam * lp.Lam + lp.Phi * lp.Phi)

	if rho != 0 {
		c := 2.0 * math.Atan2(rho, op.R2)
		sinc := math.Sin(c)
		cosc := math.Cos(c)
		lp.Lam = math.Asin(cosc * op.sinc0 + lp.Phi * sinc * op.cosc0 / rho)
		lp.Phi = math.Atan2(lp.Lam * sinc, rho * op.cosc0 * cosc - lp.Phi * op.sinc0 * sinc)
	} else {
		lp.Lam = op.phic0
		lp.Phi = 0.
	}

	// lp.Lam = op.adjust_lon(lp.Lam + op.long0)
	return lp, nil
}

// "+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 
// +y_0=463000 +ellps=bessel +towgs84=565.417,50.3319,465.552,-0.398957,0.343988,-1.8774,4.0725 
// +units=m +no_defs"
func (op *Sterea) stereaSetup(sys *core.System) error {
	PE := op.System.Ellipsoid
	// Gauss
	sys.Lam0 = sys.Lam0 * D2R
	sys.Phi0 = sys.Phi0 * D2R
	op.long0 = sys.Lam0
	op.lat0 = sys.Phi0

	op.sphi = math.Sin(op.lat0)
	op.cphi = math.Cos(op.lat0)
	op.cphi *= op.cphi
	op.rc = math.Sqrt(1.0 - PE.Es) / (1.0 - PE.Es * op.sphi * op.sphi)
	op.C = math.Sqrt(1.0 + PE.Es * op.cphi * op.cphi / (1.0 - PE.Es))
	op.phic0 = math.Asin(op.sphi / op.C)
	op.ratexp = 0.5 * op.C * PE.E
	esinp := PE.E * op.sphi
	op.K = math.Tan(0.5 * op.phic0 + FORTPI) / (
		math.Pow(math.Tan(0.5 * op.lat0 + FORTPI), op.C) * op.srat(esinp, op.ratexp))

	// Sterea
	op.sinc0 = math.Sin(op.phic0)
	op.cosc0 = math.Cos(op.phic0)
	op.R2 = 2.0 * op.rc

	return nil
}

func (op *Sterea) srat(esinp, exp float64) (float64) {
	var srat float64
	srat = math.Pow((1.0 - esinp) / (1.0 + esinp), exp)

	return srat
}

func (op *Sterea) adjust_lon(x float64) (float64) {
	if math.Abs(x) < PI {
		return x
	}

	return (x - (op.sign(x) * TWO_PI))
}

func (op *Sterea) sign(x float64) float64 {
	if x < 0.0 {
		return -1
	}

	return 1
}

func (op *Sterea) geodetic_to_geocentric(xy *core.CoordLPZ) (*core.CoordLPZ) {
	// PE := op.System.Ellipsoid
	Longitude := xy.Lam
	Latitude := xy.Phi
	// Z value not always supplied
	Height := 0.0

	/*
		* * Don't blow up if Latitude is just a little out of the value
		* * range as it may just be a rounding issue.  Also removed longitude
		* * test, it should be wrapped by cos() and sin().  NFW for PROJ.4, Sep/2001.
	*/

	if Latitude < -HALF_PI && Latitude > -1.001 * HALF_PI {
		Latitude = -HALF_PI
	} else if Latitude > HALF_PI && Latitude < 1.001 * HALF_PI {
		Latitude = HALF_PI
	} else if (Latitude < -HALF_PI) || (Latitude > HALF_PI) {
		return xy
	}

	if Longitude > PI {
		Longitude -= (2 * PI)
	}

	// sin(Latitude)
	Sin_Lat := math.Sin(Latitude)

	// cos(Latitude)
	Cos_Lat := math.Cos(Latitude)

	// Square of sin(Latitude)
	Sin2_Lat := Sin_Lat * Sin_Lat

	// Earth radius at location
	Rn := 6378137 / (math.Sqrt(1.0e0 - 0.0066943799901411 * Sin2_Lat))

	xy.Lam = (Rn + Height) * Cos_Lat * math.Cos(Longitude)
	xy.Phi = (Rn + Height) * Cos_Lat * math.Sin(Longitude)

	xy.Z = ((Rn * (1 - 0.0066943799901411)) + Height) * Sin_Lat

	return xy
}

func (op *Sterea) geocentric_from_wgs84(xy *core.CoordLPZ, datums [7]float64) (*core.CoordLPZ) {
	Dx_BF := datums[0]
	Dy_BF := datums[1]
	Dz_BF := datums[2]
	Rx_BF := datums[3]
	Ry_BF := datums[4]
	Rz_BF := datums[5]
	M_BF := datums[6]

	x_tmp := (xy.Lam - Dx_BF) / M_BF
	y_tmp := (xy.Phi - Dy_BF) / M_BF
	z_tmp := (xy.Z - Dz_BF) / M_BF

	xy.Lam = x_tmp + Rz_BF * y_tmp - Ry_BF * z_tmp
	xy.Phi = -Rz_BF * x_tmp + y_tmp + Rx_BF * z_tmp
	xy.Z = Ry_BF * x_tmp - Rx_BF * y_tmp + z_tmp

	return xy
}

func (op *Sterea) geocentric_to_geodetic(xy *core.CoordLPZ, sys *core.System) (*core.CoordLPZ) {
	PE := op.System.Ellipsoid
	genau := 1.e-12
	genau2 := (genau * genau)
	maxiter := 30

	X := xy.Lam
	Y := xy.Phi
	Z := xy.Z

	P := math.Sqrt(X * X + Y * Y)
	RR := math.Sqrt(X * X + Y * Y + Z * Z)
	Height:= 0.0
	Longitude := 0.0
	Latitude := 0.0

	// Special cases for latitude and longitude.
	if (P / 6377397.155 < genau) {
		// special case, if P=0. (X=0., Y=0.)
		// AtPole = true
		Longitude = 0.0

		// If (X,Y,Z)=(0.,0.,0.) then Height becomes semi-minor axis
		// of ellipsoid (=center of mass), Latitude becomes PI/2
		if (RR / 6377397.155 < genau) {
			Latitude = HALF_PI
			Height = PE.B * -1
			return xy
		}
	} else {
		// ellipsoidal (geodetic) longitude
		// interval: -PI < Longitude <= +PI
		Longitude = math.Atan2(Y, X)
	}

	CT := Z / RR
	ST := P / RR
	RX := 1.0 / math.Sqrt(1.0 - PE.Es * (2.0 - PE.Es) * ST * ST)
	CPHI0 := ST * (1.0 - PE.Es) * RX
	SPHI0 := CT * RX
	iter := 0
	condition := true

	CPHI := 0.0
	SPHI := 0.0

	// loop to find sin(Latitude) resp-> Latitude
	// until |sin(Latitude(iter)-Latitude(iter-1))| < genau
	for ok := true; ok; ok = condition {
		iter = iter + 1
		RN := 6377397.155 / math.Sqrt(1.0 - PE.Es * SPHI0 * SPHI0)

		/*  ellipsoidal (geodetic) height */
		Height = P * CPHI0 + Z * SPHI0 - RN * (1.0 - PE.Es * SPHI0 * SPHI0)

		RK := PE.Es * RN / (RN + Height)
		RX = 1.0 / math.Sqrt(1.0 - RK * (2.0 - RK) * ST * ST)
		CPHI = ST * (1.0 - RK) * RX
		SPHI = CT * RX
		SDPHI := SPHI * CPHI0 - CPHI * SPHI0
		CPHI0 = CPHI
		SPHI0 = SPHI
		condition = (SDPHI * SDPHI > genau2 && iter < maxiter)
	}

	// ellipsoidal (geodetic) latitude
	Latitude = math.Atan(SPHI / math.Abs(CPHI))

	xy.Lam = Longitude
	xy.Phi = Latitude
	xy.Z = Height

	return xy
}
