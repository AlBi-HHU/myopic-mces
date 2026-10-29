// Benchmark of the Java (SIRIUS structure_matching) filters, see ../README.md.
// usage: java -cp "<structure_matching>/lib/*" BenchFilters.java <pairs.csv> [values_out.csv]
import matching.algorithm.MCESDist2;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.PrintWriter;
import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

public class BenchFilters {
    public static void main(String[] args) throws Exception {
        List<String> lines = Files.readAllLines(Path.of(args[0]));
        if (lines.get(0).startsWith("smiles1,smiles2")) lines = lines.subList(1, lines.size());
        final int n = lines.size();
        final String[][] smiles = new String[n][];
        for (int i = 0; i < n; i++) {
            String[] parts = lines.get(i).split(",");
            // headerless format: id,smiles1,smiles2,...
            smiles[i] = parts.length == 2 ? parts : new String[]{parts[1], parts[2]};
        }

        SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        Method preprocess = MCESDist2.class.getDeclaredMethod("preprocessMolecule", IAtomContainer.class);
        preprocess.setAccessible(true);

        double[] f1 = new double[n], f2 = new double[n];
        // pass 0 warms up the JIT, pass 1 is reported
        for (int pass = 0; pass < 2; pass++) {
            long tParse = 0, tPre = 0, tF1 = 0, tF2 = 0;
            for (int i = 0; i < n; i++) {
                long t0 = System.nanoTime();
                IAtomContainer mol1 = parser.parseSmiles(smiles[i][0]);
                IAtomContainer mol2 = parser.parseSmiles(smiles[i][1]);
                long t1 = System.nanoTime();
                MCESDist2 dist = new MCESDist2(mol1, mol2);
                IAtomContainer p1 = (IAtomContainer) preprocess.invoke(dist, mol1);
                IAtomContainer p2 = (IAtomContainer) preprocess.invoke(dist, mol2);
                long t2 = System.nanoTime();
                f1[i] = dist.weightedDegreeBasedFilter(p1, p2);
                long t3 = System.nanoTime();
                f2[i] = dist.neighborhoodBasedFilter(p1, p2);
                long t4 = System.nanoTime();
                tParse += t1 - t0; tPre += t2 - t1; tF1 += t3 - t2; tF2 += t4 - t3;
            }
            System.out.printf("pass %d (%s): pairs %d, parse %.2fs, preprocess %.2fs%n", pass, pass == 0 ? "warmup" : "measured", n, tParse / 1e9, tPre / 1e9);
            System.out.printf("  filter1 %.2fs = %.1f us/pair%n", tF1 / 1e9, tF1 / 1e3 / n);
            System.out.printf("  filter2 %.2fs = %.1f us/pair%n", tF2 / 1e9, tF2 / 1e3 / n);
        }
        if (args.length > 1) {
            try (PrintWriter out = new PrintWriter(args[1])) {
                out.println("filter1,filter2");
                for (int i = 0; i < n; i++) out.println(f1[i] + "," + f2[i]);
            }
        }
    }
}
