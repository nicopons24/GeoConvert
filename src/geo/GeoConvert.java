package geo;

public class GeoConvert {

	/* Ellipsoid model constants (actual values here are for WGS84) */
	private static final double SEMIMAJOR = 6378137.0;
	private static final double ECCENTRICITY = 0.081819190842622;
	
	private static final int NORTH = 0;
	private static final double EAST = 500000.0;

	private static final double UMTSCALEFACTOR = 0.9996;

	/*
	 * DegToRad
	 * 
	 * Converts degrees to radians.
	 */
	public static double DegToRad(double deg) {
		return (deg / 180.0 * Math.PI);
	}

	/*
	 * RadToDeg
	 * 
	 * Converts radians to degrees.
	 */
	public static double RadToDeg(double rad) {
		return (rad / Math.PI * 180.0);
	}

	/*
	 * UTMCentralMeridian
	 * 
	 * Determines the central meridian for the given UTM zone.
	 * 
	 * Inputs: zone - An integer value designating the UTM zone, range [1,60].
	 * 
	 * Returns: The central meridian for the given UTM zone, in radians, or zero
	 * if the UTM zone parameter is outside the range [1,60]. Range of the
	 * central meridian is the radian equivalent of [-177,+177].
	 */
	public static double UTMCentralMeridian(int zone) {
		double cmeridian;

		cmeridian = DegToRad(-183.0 + (zone * 6.0));

		return cmeridian;
	}
	

	public static double[] UTMToLatLong(double north, double east, int utmZone)
    {
        // This is the lambda knot value in the reference
        double LngOrigin = UTMCentralMeridian(utmZone);

        // The following set of class constants define characteristics of the
        // ellipsoid, as defined my the WGS84 datum.  These values need to be
        // changed if a different dataum is used.

        int FalseNorth = NORTH;   // South or North?


        double Ecc = ECCENTRICITY;       // Eccentricity
        double EccSq = Ecc * Ecc;
        double Ecc2Sq = EccSq / (1. - EccSq);
        double Ecc2 = Math.sqrt(Ecc2Sq);      // Secondary eccentricity
        double E1 = ( 1 - Math.sqrt(1-EccSq) ) / ( 1 + Math.sqrt(1-EccSq) );
        double E12 = E1 * E1;
        double E13 = E12 * E1;
        double E14 = E13 * E1;

        double SemiMajor = SEMIMAJOR;         // Ellipsoidal semi-major axis (Meters)
        double FalseEast = EAST;          // UTM East bias (Meters)
        double ScaleFactor = UMTSCALEFACTOR;          // Scale at natural origin

        // Calculate the Cassini projection parameters

        double M1 = (north - FalseNorth) / ScaleFactor;
        double Mu1 = M1 / ( SemiMajor * (1 - EccSq/4.0 - 3.0*EccSq*EccSq/64.0 - 5.0*EccSq*EccSq*EccSq/256.0) );

        double Phi1 = Mu1 + (3.0*E1/2.0 - 27.0*E13/32.0) * Math.sin(2.0*Mu1)
        + (21.0*E12/16.0 - 55.0*E14/32.0)           * Math.sin(4.0*Mu1)
        + (151.0*E13/96.0)                          * Math.sin(6.0*Mu1)
        + (1097.0*E14/512.0)                        * Math.sin(8.0*Mu1);

        double sin2phi1 = Math.sin(Phi1) * Math.sin(Phi1);
        double Rho1 = (SemiMajor * (1.0-EccSq) ) / Math.pow(1.0-EccSq*sin2phi1,1.5);
        double Nu1 = SemiMajor / Math.sqrt(1.0-EccSq*sin2phi1);

        // Compute parameters as defined in the POSC specification.  T, C and D

        double T1 = Math.tan(Phi1) * Math.tan(Phi1);
        double T12 = T1 * T1;
        double C1 = Ecc2Sq * Math.cos(Phi1) * Math.cos(Phi1);
        double C12 = C1 * C1;
        double D  = (east - FalseEast) / (ScaleFactor * Nu1);
        double D2 = D * D;
        double D3 = D2 * D;
        double D4 = D3 * D;
        double D5 = D4 * D;
        double D6 = D5 * D;

        // Compute the Latitude and Longitude and convert to degrees
        double lat = Phi1 - Nu1*Math.tan(Phi1)/Rho1 * ( D2/2.0 - (5.0 + 3.0*T1 + 10.0*C1 - 4.0*C12 - 9.0*Ecc2Sq)*D4/24.0 + (61.0 + 90.0*T1 + 298.0*C1 + 45.0*T12 - 252.0*Ecc2Sq - 3.0*C12)*D6/720.0 );

        lat = RadToDeg(lat);

        double lon = LngOrigin + (D - (1.0 + 2.0*T1 + C1)*D3/6.0 + (5.0 - 2.0*C1 + 28.0*T1 - 3.0*C12 + 8.0*Ecc2Sq + 24.0*T12)*D5/120.0) / Math.cos(Phi1);

        lon = RadToDeg(lon);

        // Create a object to store the calculated Latitude and Longitude values
        double[] latLong = {lat, lon};

        // Returns a PC_LatLon object
        return latLong;
    }

}
