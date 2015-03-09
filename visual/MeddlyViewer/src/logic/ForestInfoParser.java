package logic;

import info.ForestInfo;
import info.LeafInfo;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;

import javafx.scene.chart.XYChart.Series;

import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

/**
 * A Static class which parses the JSON information from Meddly for animation.
 * 
 * @author Coleman
 *
 */
public class ForestInfoParser {
	private static ForestInfo forestInfo = null;
	private static BufferedReader br = null;
	private static JSONParser parser = null;

	/**
	 * 
	 * @throws FileNotFoundException
	 */
	public ForestInfoParser() {
		try {
			ForestInfoParser.br = new BufferedReader(new FileReader(
					"qc7.json"));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		ForestInfoParser.parser = new JSONParser();

	}

	/**
	 * 
	 * @return
	 * @throws IOException
	 */
	public static Series initalizeForestInfoFromJsonFile() throws IOException {
		try {
			String stringOfForestInfo = br.readLine();
			System.out.println("Record:\t" + stringOfForestInfo);

			Object obj;
			JSONObject jsonObject = new JSONObject();
			obj = parser.parse(stringOfForestInfo);
			jsonObject = (JSONObject) obj;
			forestInfo = new ForestInfo(jsonObject);
			return forestInfo.getSeries();

		} catch (Exception e) {
			if (br != null)
				br.close();
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * 
	 * @return
	 * @throws org.json.simple.parser.ParseException
	 * @throws IOException
	 */
	public static LeafInfo readNodeInfoFromJsonFile()
			throws org.json.simple.parser.ParseException, IOException {
		LeafInfo info = new LeafInfo();

		try {

			String sCurrentLine;
			while ((sCurrentLine = br.readLine()) != null) {
				if(br.ready() == false) break;
				System.out.println("Record:\t" + sCurrentLine);
				Object obj;
				JSONObject jsonObject = new JSONObject();
				if(!sCurrentLine.startsWith("{")) break;
				obj = parser.parse(sCurrentLine);
				jsonObject = (JSONObject) obj;
				if (jsonObject.get("f") != null) {
					Long id = (Long) jsonObject.get("f");
					System.out.println(id);
					Long level = (Long) jsonObject.get("l");
					System.out.println(level);
					Long anc = (Long) jsonObject.get("anc");
					System.out.println(anc);
					info.addAnc(anc);
					info.addLevel(level);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		br.close();
		return info;
	}

	/**
	 * 
	 * @return
	 */
	public ForestInfo getForestInfo() {
		return forestInfo;
	}
}