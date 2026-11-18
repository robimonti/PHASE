/* Read-only terrain projection for preliminary PSI quality screening.
 * Uses SNAP's own DEM/datum and orbit model. No SAR bands are resampled.
 * Output coordinates must pass an independent radar-grid check in Python.
 */
import java.io.*;
import java.lang.reflect.*;
import java.nio.file.*;
import java.util.Locale;

public final class PhaseQualityGeometry {
    private static Object geo(Class<?> type,double lat,double lon) throws Exception {
        try { return type.getConstructor(double.class,double.class).newInstance(lat,lon); }
        catch(NoSuchMethodException e) { return type.getConstructor(float.class,float.class).newInstance((float)lat,(float)lon); }
    }
    public static void main(String[] a) throws Exception {
        if(a.length!=13) throw new IllegalArgumentException("Expected source, DEM name/file/NoData, W/E/S/N, columns/rows, output, validation input/output.");
        Locale.setDefault(Locale.ROOT);
        Class<?> pc=Class.forName("org.esa.snap.core.datamodel.Product"),mc=Class.forName("org.esa.snap.core.datamodel.MetadataElement");
        Class<?> sc=Class.forName("org.jlinda.core.SLCImage"),oc=Class.forName("org.jlinda.core.Orbit"),pt=Class.forName("org.jlinda.core.Point");
        Class<?> gc=Class.forName("org.esa.snap.core.datamodel.GeoPos"),dc=Class.forName("org.esa.snap.core.dataop.dem.ElevationModel");
        Object product=Class.forName("org.esa.snap.core.dataio.ProductIO").getMethod("readProduct",File.class).invoke(null,new File(a[0]));
        if(product==null) throw new IOException("SNAP cannot read geometry metadata.");
        Object dem=null;
        try {
            Object meta=Class.forName("org.esa.snap.engine_utilities.datamodel.AbstractMetadata").getMethod("getAbstractedMetadata",pc).invoke(null,product);
            Object slc=sc.getConstructor(mc,pc).newInstance(meta,product),orbit=oc.getConstructor(mc,int.class).newInstance(meta,3);
            Method project=oc.getMethod("ell2lp",double[].class,sc),time=sc.getMethod("line2ta",double.class),satellite=oc.getMethod("getXYZ",double.class);
            Field xf=pt.getField("x"),yf=pt.getField("y"),zf=pt.getField("z");
            if(a[1].equalsIgnoreCase("External DEM")) {
                Class<?> fd=Class.forName("org.esa.snap.dem.dataio.FileElevationModel");
                try { dem=fd.getConstructor(File.class,String.class,Double.class).newInstance(new File(a[2]),"BILINEAR_INTERPOLATION",Double.valueOf(a[3])); }
                catch(NoSuchMethodException e) { dem=fd.getConstructor(File.class,String.class,double.class).newInstance(new File(a[2]),"BILINEAR_INTERPOLATION",Double.parseDouble(a[3])); }
            } else dem=Class.forName("org.esa.snap.dem.dataio.DEMFactory").getMethod("createElevationModel",String.class,String.class).invoke(null,a[1],"BILINEAR_INTERPOLATION");
            if(dem==null) throw new IOException("DEM unavailable for quality screening.");
            Method elevation=dc.getMethod("getElevation",gc);double nodata=Double.parseDouble(a[3]);
            if(!a[1].equalsIgnoreCase("External DEM")) {
                Object descriptor=dc.getMethod("getDescriptor").invoke(dem);
                nodata=((Number)Class.forName("org.esa.snap.core.dataop.dem.ElevationModelDescriptor").getMethod("getNoDataValue").invoke(descriptor)).doubleValue();
            }
            double west=Double.parseDouble(a[4]),east=Double.parseDouble(a[5]),south=Double.parseDouble(a[6]),north=Double.parseDouble(a[7]);
            int cols=Integer.parseInt(a[8]),rows=Integer.parseInt(a[9]);
            if(cols<3||rows<3||(long)cols*rows>600000||!(west<east&&south<north)) throw new IOException("Invalid screening grid.");
            try(PrintWriter w=new PrintWriter(a[10],"UTF-8")) {
                w.println("longitude,latitude,height,x,y,flat_x,flat_y,satellite_x,satellite_y,satellite_z");
                for(int row=0;row<rows;row++) {
                    double lat=south+(north-south)*row/(rows-1);
                    for(int col=0;col<cols;col++) {
                        double lon=west+(east-west)*col/(cols-1),h=((Number)elevation.invoke(dem,geo(gc,lat,lon))).doubleValue();
                        if(!Double.isFinite(h)||h==nodata) { w.printf("%.12f,%.12f,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN%n",lon,lat);continue; }
                        Object r=project.invoke(orbit,new double[]{Math.toRadians(lat),Math.toRadians(lon),h},slc);
                        Object flat=project.invoke(orbit,new double[]{Math.toRadians(lat),Math.toRadians(lon),0.0},slc);
                        Object sat=satellite.invoke(orbit,time.invoke(slc,yf.getDouble(r)));
                        w.printf("%.12f,%.12f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f%n",lon,lat,h,xf.getDouble(r),yf.getDouble(r),xf.getDouble(flat),yf.getDouble(flat),xf.getDouble(sat),yf.getDouble(sat),zf.getDouble(sat));
                    }
                }
            }
            try(BufferedReader in=Files.newBufferedReader(Paths.get(a[11]));PrintWriter w=new PrintWriter(a[12],"UTF-8")) {
                String line;while((line=in.readLine())!=null) {
                    String[] v=line.split(",");double lon=Double.parseDouble(v[0]),lat=Double.parseDouble(v[1]),h=Double.parseDouble(v[2]);
                    Object r=project.invoke(orbit,new double[]{Math.toRadians(lat),Math.toRadians(lon),h},slc);
                    w.printf("%.8f,%.8f%n",xf.getDouble(r),yf.getDouble(r));
                }
            }
        } finally {
            if(dem!=null) dc.getMethod("dispose").invoke(dem);
            pc.getMethod("dispose").invoke(product);
        }
    }
}
