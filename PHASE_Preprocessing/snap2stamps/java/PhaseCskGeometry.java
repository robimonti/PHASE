/* PHASE CSK geometry bridge. Uses the installed SNAP sensor/DEM model.
 * Reflection keeps this small Java 8 helper compatible with SNAP Java 8/21.
 * It reads metadata and elevations only; complex SAR samples are never read.
 */
import java.io.*;
import java.lang.reflect.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.Locale;

public final class PhaseCskGeometry {
    private static Object position(Class<?> type, double a, double b) throws Exception {
        try { return type.getConstructor(double.class,double.class).newInstance(a,b); }
        catch(NoSuchMethodException e) { return type.getConstructor(float.class,float.class).newInstance((float)a,(float)b); }
    }
    private static double number(Object o) { return ((Number)o).doubleValue(); }
    private static double[] axis(double low,double high,double origin,double step) throws IOException {
        if(!Double.isFinite(step)||step<=0) throw new IOException("Invalid DEM sampling interval.");
        long first=(long)Math.ceil((low-origin)/step),last=(long)Math.floor((high-origin)/step);
        long count=Math.max(0,last-first+1);
        if(count>50000000L) throw new IOException("The requested DEM grid is too large.");
        double[] values=new double[(int)count+2];values[0]=low;values[values.length-1]=high;
        for(int i=0;i<count;i++) values[i+1]=origin+(first+i)*step;
        return values;
    }
    public static void main(String[] args) throws Exception {
        if(args.length!=10) throw new IllegalArgumentException("Expected source, DEM name, DEM file, NoData, lon min/max, lat min/max, halo pixels, output JSON.");
        Locale.setDefault(Locale.ROOT);
        double west=Double.parseDouble(args[4]),east=Double.parseDouble(args[5]);
        double south=Double.parseDouble(args[6]),north=Double.parseDouble(args[7]);
        int halo=Integer.parseInt(args[8]);
        if(!(west<east && south<north && west>=-180 && east<=180 && south>-90 && north<90 && halo>=0))
            throw new IllegalArgumentException("Invalid geographic AOI or radar halo.");
        Class<?> productClass=Class.forName("org.esa.snap.core.datamodel.Product");
        Class<?> metadataClass=Class.forName("org.esa.snap.core.datamodel.MetadataElement");
        Class<?> slcClass=Class.forName("org.jlinda.core.SLCImage");
        Class<?> orbitClass=Class.forName("org.jlinda.core.Orbit");
        Class<?> geoClass=Class.forName("org.esa.snap.core.datamodel.GeoPos");
        Class<?> pointClass=Class.forName("org.jlinda.core.Point");
        Class<?> demApi=Class.forName("org.esa.snap.core.dataop.dem.ElevationModel");
        Object product=Class.forName("org.esa.snap.core.dataio.ProductIO").getMethod("readProduct",File.class).invoke(null,new File(args[0]));
        if(product==null) throw new IOException("SNAP cannot read the CSK source product: "+args[0]);
        Object dem=null;
        try {
            Object metadata=Class.forName("org.esa.snap.engine_utilities.datamodel.AbstractMetadata")
                    .getMethod("getAbstractedMetadata",productClass).invoke(null,product);
            Object slc=slcClass.getConstructor(metadataClass,productClass).newInstance(metadata,product);
            Object orbit=orbitClass.getConstructor(metadataClass,int.class).newInstance(metadata,3);
            int width=((Number)productClass.getMethod("getSceneRasterWidth").invoke(product)).intValue();
            int height=((Number)productClass.getMethod("getSceneRasterHeight").invoke(product)).intValue();
            if(args[1].equalsIgnoreCase("External DEM")) {
                Class<?> fileDem=Class.forName("org.esa.snap.dem.dataio.FileElevationModel");
                try { dem=fileDem.getConstructor(File.class,String.class,Double.class)
                        .newInstance(new File(args[2]),"BILINEAR_INTERPOLATION",Double.valueOf(args[3])); }
                catch(NoSuchMethodException e) { dem=fileDem.getConstructor(File.class,String.class,double.class)
                        .newInstance(new File(args[2]),"BILINEAR_INTERPOLATION",Double.parseDouble(args[3])); }
                // Keep SNAP's default vertical-datum handling, as in TopoPhaseRemoval.
            } else {
                dem=Class.forName("org.esa.snap.dem.dataio.DEMFactory")
                        .getMethod("createElevationModel",String.class,String.class).invoke(null,args[1],"BILINEAR_INTERPOLATION");
            }
            if(dem==null) throw new IOException("SNAP cannot initialize the selected DEM.");
            Class<?> pixelClass=Class.forName("org.esa.snap.core.datamodel.PixelPos");
            Method geoPosition=demApi.getMethod("getGeoPos",pixelClass);
            Object origin=geoPosition.invoke(dem,position(pixelClass,.5,.5));
            Object eastPixel=geoPosition.invoke(dem,position(pixelClass,1.5,.5));
            Object southPixel=geoPosition.invoke(dem,position(pixelClass,.5,1.5));
            Method getLat=geoClass.getMethod("getLat"),getLon=geoClass.getMethod("getLon");
            double originLat=number(getLat.invoke(origin)),originLon=number(getLon.invoke(origin));
            double dx=Math.abs(number(getLon.invoke(eastPixel))-originLon);
            double dy=Math.abs(number(getLat.invoke(southPixel))-originLat);
            if(Math.abs(number(getLat.invoke(eastPixel))-originLat)>1e-9 ||
                    Math.abs(number(getLon.invoke(southPixel))-originLon)>1e-9)
                throw new IOException("The external DEM must use a north-up geographic WGS84 grid.");
            // Include actual DEM pixel centres and all AOI edges, not just corners.
            double[] longitudes=axis(west,east,originLon,dx),latitudes=axis(south,north,originLat,dy);
            int nx=longitudes.length-1,ny=latitudes.length-1;
            long total=(long)longitudes.length*latitudes.length;
            if(total>50000000L) throw new IOException("AOI/DEM geometry scan exceeds 50 million samples; reduce the AOI or provide an appropriately sampled DEM.");
            Method elevation=demApi.getMethod("getElevation",geoClass);
            Method inverse=orbitClass.getMethod("ell2lp",double[].class,slcClass);
            Field xField=pointClass.getField("x"),yField=pointClass.getField("y");
            double noData=Double.parseDouble(args[3]);
            if(!args[1].equalsIgnoreCase("External DEM")) {
                Object descriptor=demApi.getMethod("getDescriptor").invoke(dem);
                noData=number(Class.forName("org.esa.snap.core.dataop.dem.ElevationModelDescriptor")
                        .getMethod("getNoDataValue").invoke(descriptor));
            }
            double minX=Double.POSITIVE_INFINITY,minY=minX,maxX=Double.NEGATIVE_INFINITY,maxY=maxX;
            double minH=Double.POSITIVE_INFINITY,maxH=Double.NEGATIVE_INFINITY;
            long missing=0,invalid=0,inside=0;
            for(int row=0;row<=ny;row++) {
                double lat=latitudes[row];
                for(int col=0;col<=nx;col++) {
                    double lon=longitudes[col];
                    double z=number(elevation.invoke(dem,position(geoClass,lat,lon)));
                    if(!Double.isFinite(z)||z==noData) { missing++;continue; }
                    Object radar=inverse.invoke(orbit,new double[]{Math.toRadians(lat),Math.toRadians(lon),z},slc);
                    // ell2lp returns zero-based sample/line coordinates used by SNAP's radar coding.
                    double x=xField.getDouble(radar),y=yField.getDouble(radar);
                    if(!Double.isFinite(x)||!Double.isFinite(y)) { invalid++;continue; }
                    minX=Math.min(minX,x);maxX=Math.max(maxX,x);minY=Math.min(minY,y);maxY=Math.max(maxY,y);
                    minH=Math.min(minH,z);maxH=Math.max(maxH,z);
                    if(x>=0&&x<width&&y>=0&&y<height) inside++;
                }
                if(row==0||row==ny||row%Math.max(1,ny/10)==0)
                    System.out.println("CSK terrain-aware AOI: "+Math.round(100.0*row/ny)+"% (DEM and orbit).");
            }
            if(missing>0) throw new IOException("The DEM has "+missing+" missing elevations among "+total+" AOI samples. Provide complete elevation coverage; missing heights are not replaced with zero.");
            if(invalid>0) throw new IOException("The SNAP orbit model failed for "+invalid+" AOI samples.");
            if(inside==0) throw new IOException("The requested ground AOI does not intersect this CSK acquisition after terrain-aware projection.");
            int x0=Math.max(0,(int)Math.floor(minX)-halo),y0=Math.max(0,(int)Math.floor(minY)-halo);
            int x1=Math.min(width-1,(int)Math.ceil(maxX)+halo),y1=Math.min(height-1,(int)Math.ceil(maxY)+halo);
            if(x1<x0||y1<y0) throw new IOException("Empty terrain-aware radar crop.");
            String json=String.format(Locale.ROOT,
                "{\"schema\":1,\"algorithm\":\"terrain-aware-v1\",\"requestedBounds\":[%.12f,%.12f,%.12f,%.12f],\"region\":[%d,%d,%d,%d],\"sourceSize\":[%d,%d],\"haloPixels\":%d,\"demSamples\":%d,\"insideSourceSamples\":%d,\"heightRange\":[%.8f,%.8f],\"unclippedRadarBounds\":[%.8f,%.8f,%.8f,%.8f]}",
                west,east,south,north,x0,y0,x1-x0+1,y1-y0+1,width,height,halo,total,inside,minH,maxH,minX,maxX,minY,maxY);
            Files.write(Paths.get(args[9]),json.getBytes(StandardCharsets.UTF_8));
            if(inside<total) System.out.println("Requested AOI extends beyond the CSK scene; processing its available intersection ("+inside+"/"+total+" ground samples).");
        } finally {
            if(dem!=null) demApi.getMethod("dispose").invoke(dem);
            productClass.getMethod("dispose").invoke(product);
        }
    }
}
