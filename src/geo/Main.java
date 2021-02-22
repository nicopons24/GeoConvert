package geo;

public class Main {

	public static void main(String[] args) {
		
		double[] latlong = new double[2];
		double x = 4372961.315235150977969; //Point[1]
		double y = 725514.652338833780959;//Point[0]
		int zone = 30;
		
		latlong = GeoConvert.UTMToLatLong(x, y, zone);
		System.out.println(latlong[0]);
		System.out.println(latlong[1]);
		
	}

}
