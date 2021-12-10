package proteomics;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Library.Binomial;
import proteomics.FDR.EstimateFDR;
import proteomics.Segment.InferenceSegment;
import proteomics.Spectrum.FilterSpectra;
import proteomics.Types.*;
import proteomics.PTM.FindPTM;
import proteomics.Index.BuildIndex;
import proteomics.Parameter.Parameter;
import proteomics.Spectrum.PreSpectra;
import proteomics.Library.PrepareSpectrum;
import proteomics.Library.MassTool;
//import proteomics.TheoSeq.MassTool;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLParsingException;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.sql.*;
import java.text.SimpleDateFormat;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantLock;

public class PIPI {

    private static final Logger logger = LoggerFactory.getLogger(PIPI.class);
    private static final double C13Diff = 1.00335483f;
    public static final int[] debugScanNumArray = new int[]{};
    public static final String versionStr = "1.5.0";
    static final boolean useXcorr = true;

    public static final boolean DEV = false;

    public static void main(String args[]) {
        // Process inputs
        long startTime = System.nanoTime();
        if (args.length != 2) {
            help();
        }

        // Set parameters
        String parameterPath = args[0].trim();
        String spectraPath = args[1].trim();

        logger.info("Running PIPI version {}.", versionStr);

        String dbName = null;
        String hostName = "unknown-host";

        try {
            hostName = InetAddress.getLocalHost().getHostName();
            logger.info("Computer: {}.", hostName);
        } catch (UnknownHostException ex) {
            logger.warn("Cannot get the computer's name.");
        }

        try {
            logger.info("Spectra: {}, parameter: {}.", spectraPath, parameterPath);

            dbName = String.format(Locale.US, "PIPI.%s.%s.temp.db", hostName, new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss").format(Calendar.getInstance().getTime()));
            new PIPI(parameterPath, spectraPath, dbName);
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
        } finally {
            if (dbName != null) {
                (new File(dbName)).delete();
                (new File(dbName + "-wal")).delete();
                (new File(dbName + "-shm")).delete();
            }
        }

        double totalHour = (double) (System.nanoTime() - startTime) * 1e-9 / 3600;
        logger.info("Running time: {} hours.", totalHour);
        logger.info("Done!");
    }

    private PIPI(String parameterPath, String spectraPath, String dbName) throws Exception {
        // Get the parameter map
        Parameter parameter = new Parameter(parameterPath);
        Map<String, String> parameterMap = parameter.returnParameterMap();
        double ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        double ms1Tolerance = Double.valueOf(parameterMap.get("ms1_tolerance"));
        double leftInverseMs1Tolerance = 1 / (1 + ms1Tolerance * 1e-6);
        double rightInverseMs1Tolerance = 1 / (1 - ms1Tolerance * 1e-6);
        double minClear = Double.valueOf(parameterMap.get("min_clear_mz"));
        double maxClear = Double.valueOf(parameterMap.get("max_clear_mz"));
        int ms1ToleranceUnit = Integer.valueOf(parameterMap.get("ms1_tolerance_unit"));
        int maxMs2Charge = Integer.valueOf(parameterMap.get("max_ms2_charge"));
        int batchSize = Integer.valueOf(parameterMap.get("batch_size"));
        String percolatorPath = parameterMap.get("percolator_path");
        boolean outputPercolatorInput = (Integer.valueOf(parameterMap.get("output_percolator_input")) == 1);

        // print all the parameters
        logger.info("Parameters:");
        for (String k : parameterMap.keySet()) {
            logger.info("{} = {}", k, parameterMap.get(k));
        }

        // Check if Percolator can be executed.
        if (!(new File(percolatorPath)).exists()) {
            throw new NullPointerException(String.format(Locale.US, "Cannot find Percolator from %s.", percolatorPath));
        }
        if (!(new File(percolatorPath)).canExecute()) {
            throw new Exception(String.format(Locale.US, "Percolator (%s) exits but cannot be executed.", percolatorPath));
        }

        String[] tempArray = parameterMap.get("ms_level").split(",");
        Set<Integer> msLevelSet = new HashSet<>(tempArray.length + 1, 1);
        for (String temp : tempArray) {
            msLevelSet.add(Integer.valueOf(temp));
        }

        String labelling = "N14";
        if (parameterMap.get("15N").trim().contentEquals("1")) {
            logger.info("N15 mode is on...");
            labelling = "N15";
        }

        if (parameterMap.get("add_decoy").contentEquals("0")) {
            logger.warn("add_decoy = 0. Won't search the decoy sequences and estimate FDR.");
        }

        if (parameterMap.get("add_contaminant").contentEquals("0")) {
            logger.warn("add_contaminant = 0. Won't search the build-in contaminant proteins.");
        }

        logger.info("Indexing protein database...");
        BuildIndex buildIndexObj = new BuildIndex(parameterMap, labelling);
        MassTool massToolObj = buildIndexObj.returnMassToolObj();

        logger.info("Reading PTM database...");
//        Map<String, TreeSet<Integer>> siteMass1000Map = readPTMDb(parameterMap, buildIndexObj.getFixModMap());

//        int tempMax = 0;
//        int tempMin = 99999;
//        for (TreeSet<Integer> tempSet : siteMass1000Map.values()) {
//            if (tempSet.last() > tempMax) {
//                tempMax = tempSet.last();
//            }
//            if (tempSet.first() < tempMin) {
//                tempMin = tempSet.first();
//            }
//        }

        double minPtmMass = Double.valueOf(parameterMap.get("min_ptm_mass"));
        double maxPtmMass = Double.valueOf(parameterMap.get("max_ptm_mass"));

        logger.info("Reading spectra...");
        File spectraFile = new File(spectraPath);
        if ((!spectraFile.exists() || (spectraFile.isDirectory()))) {
            throw new FileNotFoundException("The spectra file not found.");
        }
        String[] temp = spectraPath.split("\\.");
        String ext = temp[temp.length - 1];
        JMzReader spectraParser;
        if (ext.contentEquals("mzXML")) {
            spectraParser = new MzXMLFile(spectraFile);
        } else if (ext.toLowerCase().contentEquals("mgf")) {
            spectraParser = new MgfFile(spectraFile);
        } else {
            throw new Exception(String.format(Locale.US, "Unsupported file format %s. Currently, PIPI only support mzXML and MGF.", ext));
        }

        String sqlPath = "jdbc:sqlite:" + dbName;
        Class.forName("org.sqlite.JDBC").newInstance();

        FilterSpectra filterSpectraObj = new FilterSpectra(spectraParser, parameterMap, massToolObj.getMassTable());
        FilterSpectra.MassScan[] scanNumArray = filterSpectraObj.getScanNumArray();

        logger.info("Useful MS/MS spectra number: {}", scanNumArray.length);
        logger.info("Start searching...");
        int threadNum = Integer.valueOf(parameterMap.get("thread_num"));
        if (threadNum == 0) {
            threadNum = 3 + Runtime.getRuntime().availableProcessors();
        }
        if (debugScanNumArray.length > 0) {
            threadNum = 1;
        }
        ExecutorService threadPool = Executors.newFixedThreadPool(threadNum);

        PrepareSpectrum preSpectrum = new PrepareSpectrum(massToolObj);
//        PreSpectra preSpectraObj = new PreSpectra(spectraParser, parameterMap, massToolObj, ext, scanNumArray, 0, scanNumArray.length);
//        Map<Integer, SpectrumEntry> numSpectrumMap = preSpectraObj.returnNumSpectrumMap();
//        TreeMap<Double, List<Integer>> massNumMap = preSpectraObj.returnMassNumMap();
        ArrayList<Future<PIPIWrap.Entry>> taskList = new ArrayList<>(scanNumArray.length + 10);
        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanId, precursorCharge, precursorMass FROM spectraTable");
        ReentrantLock lock = new ReentrantLock();
        Binomial binomial = new Binomial(Integer.valueOf(parameterMap.get("max_peptide_length")) * 2);
        InferenceSegment inference3SegmentObj = new InferenceSegment(buildIndexObj, ms2Tolerance, 3);
        while (sqlResultSet.next()) {
            String scanId = sqlResultSet.getString("scanId");
            int precursorCharge = sqlResultSet.getInt("precursorCharge");
            double precursorMass = sqlResultSet.getDouble("precursorMass");
            taskList.add(threadPool.submit(new PIPIWrap(inference3SegmentObj, buildIndexObj, massToolObj, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance, ms1ToleranceUnit, ms2Tolerance, minPtmMass, maxPtmMass, Math.min(precursorCharge > 1 ? precursorCharge - 1 : 1, 3), spectraParser, minClear, maxClear, lock, scanId, precursorCharge, precursorMass, preSpectrum, sqlPath, binomial)));


        }










//        InferenceSegment inference3SegmentObj = new InferenceSegment(buildIndexObj, ms2Tolerance, 3);

        List<FinalResultEntry> finalScoredPsms = new LinkedList<>();
        Map<Integer, ChargeMassTuple> numChargeMassMap = new HashMap<>();
        int startIdx;
        int endIdx = 0;
        while (endIdx < scanNumArray.length) {
            startIdx = endIdx;
            if (batchSize == 0) {
                endIdx = scanNumArray.length;
            } else {
                endIdx = Math.min(startIdx + batchSize, scanNumArray.length);
            }

//            logger.info("Searching batch {} - {} ({}%)...", startIdx, endIdx, String.format("%.1f", (double) endIdx * 100 / (double) scanNumArray.length));
            PreSpectra preSpectraObj = new PreSpectra(spectraParser, parameterMap, massToolObj, ext, scanNumArray, startIdx, endIdx);
            Map<Integer, SpectrumEntry> numSpectrumMap = preSpectraObj.returnNumSpectrumMap();
            TreeMap<Double, List<Integer>> massNumMap = preSpectraObj.returnMassNumMap();
            numChargeMassMap.putAll(preSpectraObj.getNumChargeMassMap());

            if (!massNumMap.isEmpty()) {
                int batchSize2 = (massNumMap.size() / threadNum) + 1;
                Double[] massArray = massNumMap.keySet().toArray(new Double[massNumMap.size()]);
//                Collection<PIPIWrap> taskList = new LinkedList<>();
                for (int i = 0; i < threadNum; ++i) {
                    int leftIdx = i * batchSize2;
                    int rightIdx = Math.min((i + 1) * batchSize2, massArray.length - 1);
                    if (leftIdx > massArray.length - 1) {
                        break;
                    }

                    if (rightIdx < massArray.length - 1) {
                        NavigableMap<Double, List<Integer>> subMassNumMap = massNumMap.subMap(massArray[leftIdx], true, massArray[rightIdx], false);
//                        taskList.add(new PIPIWrap(buildIndexObj, massToolObj, inference3SegmentObj, numSpectrumMap, subMassNumMap, siteMass1000Map, ms1Tolerance, ms1ToleranceUnit, ms2Tolerance, minPtmMass, maxPtmMass, maxMs2Charge, startIdx));
                    } else {
                        NavigableMap<Double, List<Integer>> subMassNumMap = massNumMap.subMap(massArray[leftIdx], true, massArray[rightIdx], true);
//                        taskList.add(new PIPIWrap(buildIndexObj, massToolObj, inference3SegmentObj, numSpectrumMap, subMassNumMap, siteMass1000Map, ms1Tolerance, ms1ToleranceUnit, ms2Tolerance, minPtmMass, maxPtmMass, maxMs2Charge, startIdx));
                    }
                }
                try {
                    List<Future<List<FinalResultEntry>>> tempResultList = threadPool.invokeAll(taskList);
                    for (Future<List<FinalResultEntry>> tempResult : tempResultList) {
                        if (tempResult.isDone() && !tempResult.isCancelled()) {
                            finalScoredPsms.addAll(tempResult.get());
                        } else {
                            logger.error("Threads were not finished normally.");
                            System.exit(1);
                        }
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                    logger.error(ex.getMessage());
                    System.exit(1);
                }
            }
        }

        // shutdown threads.
        threadPool.shutdown();
        try {
            if (!threadPool.awaitTermination(60, TimeUnit.SECONDS)) {
                threadPool.shutdownNow();
                if (!threadPool.awaitTermination(60, TimeUnit.SECONDS))
                    System.err.println("Pool did not terminate");
            }
        } catch (InterruptedException ie) {
            threadPool.shutdownNow();
            Thread.currentThread().interrupt();
            logger.error("Threads were not finished normally.");
            System.exit(1);
        }

        logger.info("Estimating FDR...");
        // estimate T-D FDR
        new EstimateFDR(finalScoredPsms);

        // estimate Percolator FDR
        String percolatorInputFileName = spectraPath + ".input.temp";
        String percolatorOutputFileName = spectraPath + ".output.temp";
        Map<String, Set<String>> peptideProteinMap = buildIndexObj.returnPepProMap();
        Map<String, String> decoyPeptideProteinMap = buildIndexObj.returnDecoyPepProMap();
        writePercolator(finalScoredPsms, numChargeMassMap, peptideProteinMap, decoyPeptideProteinMap, percolatorInputFileName, buildIndexObj.returnFixModMap());
        Map<Integer, PercolatorEntry> percolatorResultMap = runPercolator(percolatorPath, percolatorInputFileName, percolatorOutputFileName);

        if (percolatorResultMap.isEmpty()) {
            logger.info("Percolator failed to estimate FDR. The results won't contain percolator_score, posterior_error_prob, and percolator_q_value.");
        }

        if (!outputPercolatorInput) {
            (new File(percolatorInputFileName)).delete();
            (new File(percolatorOutputFileName)).delete();
        }

        logger.info("Saving results...");
        writeFinalResult(finalScoredPsms, percolatorResultMap, numChargeMassMap, peptideProteinMap, spectraPath + ".pipi.csv", buildIndexObj.returnFixModMap());

        logger.info("Done.");
    }

    private static void help() {
        String helpStr = "PIPI version " + versionStr + "\r\n"
                + "A tool identifying peptides with unlimited PTM.\r\n"
                + "Author: Fengchao Yu\r\n"
                + "Email: fyuab@connect.ust.hk\r\n"
                + "ECL usage: java -Xmx25g -jar /path/to/PIPI.jar <parameter_file> <data_file>\r\n"
                + "\t<parameter_file>: parameter file. Can be download along with PIPI.\r\n"
                + "\t<data_file>: spectra data file (mzXML)\r\n"
                + "\texample: java -Xmx32g -jar ECL.jar parameter.def data.mzxml\r\n";
        System.out.print(helpStr);
        System.exit(1);
    }

    private static void writePercolator(List<FinalResultEntry> finalScoredResult, Map<Integer, ChargeMassTuple> numChargeMassMap, Map<String, Set<String>> peptideProteinMap, Map<String, String> decoyPeptideProteinMap, String resultPath, Map<String, Double> fixModMap) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(resultPath))) {
            writer.write("id\tlabel\tscannr\txcorr\tdelta_c\tdelta_L_c\tnegative_log10_e_value\tnormalized_cross_corr\tglobal_search_rank\tabs_ppm\tIonFrac\tmatched_high_peak_frac\tcharge1\tcharge2\tcharge3\tcharge4\tcharge5\tcharge6\tvar_PTM_num\tpeptide\tprotein\n");
            for (FinalResultEntry entry : finalScoredResult) {
                double expMass = numChargeMassMap.get(entry.getScanNum()).mass;
                Peptide peptide = entry.getPeptide();
                double theoMass = peptide.getPrecursorMass();
                double massDiff = getMassDiff(expMass, theoMass, C13Diff);
                String proteinIdStr = "";

                if (!entry.isDecoy()) {
                    for (String temp : peptideProteinMap.get(peptide.getPTMFreeSeq())) {
                        proteinIdStr += temp + ";";
                    }
                } else {
                    proteinIdStr = decoyPeptideProteinMap.get(peptide.getPTMFreeSeq());
                }

                StringBuilder sb = new StringBuilder(20);
                int charge = numChargeMassMap.get(entry.getScanNum()).charge;
                for (int i = 0; i < 6; ++i) {
                    if (i == charge - 1) {
                        sb.append(1);
                    } else {
                        sb.append(0);
                    }
                    sb.append("\t");
                }

                if (entry.isDecoy()) {
                    writer.write(entry.getScanNum() + "\t-1\t" + entry.getScanNum() + "\t" + entry.getScore() + "\t" + entry.getDeltaC() + "\t" + entry.getDeltaLC() + "\t" + entry.getNegativeLog10EValue() + "\t" + entry.getNormalizedCrossXcorr() + "\t" + entry.getGlobalSearchRank() + "\t" + Math.abs(massDiff * 1e6f / theoMass) + "\t" + entry.getIonFrac() + "\t" + entry.getMatchedHighestIntensityFrac() + "\t" + sb.toString() + peptide.getVarPTMNum() + "\t" + peptide.getLeftFlank() + "." + peptide.getPTMContainedString(fixModMap) + "." + peptide.getRightFlank() + "\t" + proteinIdStr + "\n");
                } else {
                    writer.write(entry.getScanNum() + "\t1\t" + entry.getScanNum() + "\t" + entry.getScore() + "\t" + entry.getDeltaC() + "\t" + entry.getDeltaLC() + "\t" + entry.getNegativeLog10EValue() + "\t" + entry.getNormalizedCrossXcorr() + "\t" + entry.getGlobalSearchRank() + "\t" + Math.abs(massDiff * 1e6f / theoMass) + "\t" + entry.getIonFrac() + "\t" + entry.getMatchedHighestIntensityFrac() + "\t" + sb.toString() + peptide.getVarPTMNum() + "\t" + peptide.getLeftFlank() + "." + peptide.getPTMContainedString(fixModMap) + "." + peptide.getRightFlank() + "\t" + proteinIdStr + "\n");
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    private static Map<Integer, PercolatorEntry> runPercolator(String percolatorPath, String percolatorInputFileName, String percolatorOutputFileName) {
        Map<Integer, PercolatorEntry> percolatorResultMap = new HashMap<>();
        try {
            if ((new File(percolatorPath)).exists()) {
                Process ps = Runtime.getRuntime().exec(percolatorPath + " --only-psms --results-psms " + percolatorOutputFileName + " " + percolatorInputFileName);
                ps.waitFor();

                if (!(new File(percolatorOutputFileName).exists())) {
                    logger.warn("Error in running Percolator. The results won't contain percolator_score, posterior_error_prob, and percolator_q_value.");
                    return percolatorResultMap;
                }

                BufferedReader reader = new BufferedReader(new FileReader(percolatorOutputFileName));
                String line;
                while ((line = reader.readLine()) != null) {
                    line = line.trim();
                    if (!line.startsWith("PSMId")) {
                        String[] parts = line.split("\t");
                        percolatorResultMap.put(Integer.valueOf(parts[0]), new PercolatorEntry(Double.valueOf(parts[1]), parts[2], parts[3]));
                    }
                }
                reader.close();

                if (DEV) {
                    reader = new BufferedReader(new InputStreamReader(ps.getErrorStream()));
                    while ((line = reader.readLine()) != null) {
                        System.err.print(line);
                    }
                    reader.close();
                }
            } else {
                logger.error("Cannot find Percolator for estimating Percolator Q-Value. The results won't contain percolator_score, posterior_error_prob, and percolator_q_value.");
                return percolatorResultMap;
            }
        } catch (IOException | InterruptedException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            return percolatorResultMap;
        }

        return percolatorResultMap;
    }

    private static void writeFinalResult(List<FinalResultEntry> finalScoredPsms, Map<Integer, PercolatorEntry> percolatorResultMap, Map<Integer, ChargeMassTuple> numChainMassMap, Map<String, Set<String>> peptideProteinMap, String outputPath, Map<String, Double> fixModMap) {
        TreeMap<Double, List<String>> tempMap = new TreeMap<>();
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath))) {
            writer.write("scan_num,peptide,charge,theo_mass,exp_mass,ppm,protein_ID,xcorr,e_value,naive_q_value,percolator_score,posterior_error_prob,percolator_q_value\n");
            for (FinalResultEntry entry : finalScoredPsms) {
                if (!entry.isDecoy()) {
                    int scanNum = entry.getScanNum();
                    double expMass = numChainMassMap.get(scanNum).mass;
                    int charge = numChainMassMap.get(scanNum).charge;
                    Peptide peptide = entry.getPeptide();
//                    System.out.println(scanNum+" "+peptide.getPTMFreeSeq());
                    double theoMass = peptide.getPrecursorMass();
                    double massDiff = getMassDiff(expMass, theoMass, C13Diff);
                    double ppm = Math.abs(massDiff * 1e6f / theoMass);


                    String proteinIdStr = "";
                    for (String tempStr : peptideProteinMap.get(peptide.getPTMFreeSeq())) {
                        proteinIdStr += tempStr + ";";
                    }

                    String str;
                    boolean sortedByPercolatorScore = true;
                    if (percolatorResultMap.containsKey(scanNum)) {
                        PercolatorEntry percolatorEntry = percolatorResultMap.get(scanNum);
                        str = String.format("%d,%s,%d,%.4f,%.4f,%.2f,%s,%.4f,%E,%.4f,%.3f,%s,%s\n", scanNum, peptide.getPTMContainedString(fixModMap), charge, theoMass, expMass, ppm, proteinIdStr, entry.getScore(), entry.getEValue(), entry.getQValue(), percolatorEntry.percolatorScore, percolatorEntry.PEP, percolatorEntry.qValue);
                    } else {
                        sortedByPercolatorScore = false;
                        str = String.format("%d,%s,%d,%.4f,%.4f,%.2f,%s,%.4f,%E,%.4f,%s,%s,%s\n", scanNum, peptide.getPTMContainedString(fixModMap), charge, theoMass, expMass, ppm, proteinIdStr, entry.getScore(), entry.getEValue(), entry.getQValue() , "-", "-", "-");
                    }

                    if (sortedByPercolatorScore) {
                        if (tempMap.containsKey(percolatorResultMap.get(scanNum).percolatorScore)) {
                            tempMap.get(percolatorResultMap.get(scanNum).percolatorScore).add(str);
                        } else {
                            List<String> tempList = new LinkedList<>();
                            tempList.add(str);
                            tempMap.put(percolatorResultMap.get(scanNum).percolatorScore, tempList);
                        }
                    } else {
                        if (tempMap.containsKey(entry.getNegativeLog10EValue())) {
                            tempMap.get(entry.getNegativeLog10EValue()).add(str);
                        } else {
                            List<String> tempList = new LinkedList<>();
                            tempList.add(str);
                            tempMap.put(entry.getNegativeLog10EValue(), tempList);
                        }
                    }
                }
            }

            Double[] tempArray = tempMap.keySet().toArray(new Double[tempMap.size()]);
            for (int i = tempArray.length - 1; i >= 0; --i) {
                List<String> tempList = tempMap.get(tempArray[i]);
                for (String tempStr : tempList) {
                    writer.write(tempStr);
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    private static double getMassDiff(double expMass, double theoMass, double C13Diff) {
        double massDiff1 = expMass - theoMass;
        double massDiff2 = expMass - theoMass - C13Diff;
        double massDiff3 = expMass - theoMass - 2 * C13Diff;
        double absMassDiff1 = Math.abs(massDiff1);
        double absMassDiff2 = Math.abs(massDiff2);
        double absMassDiff3 = Math.abs(massDiff3);

        if ((absMassDiff1 <= absMassDiff2) && (absMassDiff1 <= absMassDiff2)) {
            return massDiff1;
        } else if ((absMassDiff2 <= absMassDiff1) && (absMassDiff2 <= absMassDiff3)) {
            return massDiff2;
        } else {
            return massDiff3;
        }
    }

    private Map<String, TreeSet<Integer>> readPTMDb(Map<String, String> parameterMap, Map<String, Double> fixModMap) {
        Map<String, TreeSet<Integer>> siteMass1000Map = new HashMap<>();
        double minPtmMass = Double.valueOf(parameterMap.get("min_ptm_mass"));
        double maxPtmMass = Double.valueOf(parameterMap.get("max_ptm_mass"));

        try (BufferedReader reader = new BufferedReader(new FileReader(parameterMap.get("PTM_db")))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.startsWith("site") || line.startsWith("#")) {
                    continue;
                }

                String[] parts = line.split("\t");
                String site = parts[0];
                String position = parts[3];
                double mass = Double.valueOf(parts[2]);

                if ((mass > maxPtmMass) || mass < minPtmMass) {
                    continue;
                }

                int mass1000 = (int) Math.floor(mass * 1000);

                String siteString;
                if (site.contentEquals("N-term") && (Math.abs(fixModMap.get("n")) < 1e-6)) {
                    if (position.contentEquals("PROTEIN_N")) {
                        siteString = "PROTEIN_N";
                    } else {
                        siteString = "PEPTIDE_N";
                    }
                } else if (site.contentEquals("C-term")) {
                    if (position.contentEquals("PROTEIN_C")) {
                        siteString = "PROTEIN_C";
                    } else {
                        siteString = "PEPTIDE_C";
                    }
                } else if (Math.abs(fixModMap.get(site)) < 1e-6) { // fix modified amino acid cannot be modified again.
                    if (position.contentEquals("PROTEIN_N")) {
                        siteString = site + "-PROTEIN_N";
                    } else if (position.contentEquals("PROTEIN_C")) {
                        siteString = site + "-PROTEIN_C";
                    } else if (position.contentEquals("PEPTIDE_N")) {
                        siteString = site + "-PEPTIDE_N";
                    } else if (position.contentEquals("PEPTIDE_C")) {
                        siteString = site + "-PEPTIDE_C";
                    } else {
                        siteString = site;
                    }
                } else {
                    continue;
                }

                if (siteMass1000Map.containsKey(siteString)) {
                    siteMass1000Map.get(siteString).add(mass1000);
                } else {
                    TreeSet<Integer> tempSet = new TreeSet<>();
                    tempSet.add(mass1000);
                    siteMass1000Map.put(siteString, tempSet);
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }

        return siteMass1000Map;
    }
}
