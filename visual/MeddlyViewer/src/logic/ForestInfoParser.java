
// $Id$

/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2009, Iowa State University Research Foundation, Inc.

    This library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published 
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this library.  If not, see <http://www.gnu.org/licenses/>.
*/

package logic;

import info.ForestInfo;
import info.LeafInfo;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;

import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.chart.BarChart;
import javafx.scene.chart.XYChart;
import javafx.scene.chart.XYChart.Series;
import javafx.util.Duration;

import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

/**
 * A Static class which parses the information from Meddly for animation.
 * 
 * @author Coleman
 *
 */
public class ForestInfoParser {
	private static ForestInfo forestInfo = null;
	private static BufferedReader br = null;

	/**
	 * 
	 * @throws FileNotFoundException
	 */
	public ForestInfoParser() {
		try {
			ForestInfoParser.br = new BufferedReader(new FileReader("qc3.txt"));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * 
	 * @return
	 * @throws IOException
	 * 
	 */
	public static Series initalizeForestInfo() throws IOException {
		try {

			// FILETYPE PARSE BLOCK
			String fileType = br.readLine(); // does nothing yet

			br.readLine(); // consume known comment line

			System.out.println(fileType);// error checking

			String stringOfForestInfo = br.readLine();
			System.out.println("Record: " + stringOfForestInfo);
			// END OF FILETYPE PARSE BLOCK

			// ID PARSE BLOCK
			int id = (int) stringOfForestInfo.charAt(2);
			// END OF ID PARSE BLOCK

			// NAME PARSE BLOCK
			String name = "";
			int index = -1;
			if (stringOfForestInfo.charAt(4) == '"') {
				index = 5;
				while (stringOfForestInfo.charAt(index) != '"') {
					name += stringOfForestInfo.charAt(index);
					index++;
				}
				// If there isn't a " at index 4, then there is no name in the
				// file, so name it untitled, and continue.
			} else {
				name = "Untitled";
				index = 4; // index used to get locations of subsequent items
							// after the name of the tree.
			}
			// END OF NAME PARSE BLOCK

			// LEFT AND RIGHT INITIAL NODE PARSE BLOCK
			if (index != 4)
				index += 2; // if the index is not 4, then name was parsed and
							// we are on a ", so increment to the next char
							// index.
			// Use a string builder to build the strings for left and right
			// counts, length of string unknown.
			StringBuilder leftAndRightCountStringConverter = new StringBuilder();
			leftAndRightCountStringConverter.append(stringOfForestInfo
					.charAt(index));
			while (stringOfForestInfo.charAt(index + 1) != ' ') {
				index++;
				leftAndRightCountStringConverter.append(stringOfForestInfo
						.charAt(index));
			}
			int leftCount = Integer.parseInt(leftAndRightCountStringConverter
					.toString());
			leftAndRightCountStringConverter.delete(0,
					leftAndRightCountStringConverter.length());

			// End of left node parsing, beginning of right node parsing
			leftAndRightCountStringConverter.append(stringOfForestInfo
					.charAt(index += 2));
			while (stringOfForestInfo.charAt(index + 1) != ' ') {
				index++;
				leftAndRightCountStringConverter.append(stringOfForestInfo
						.charAt(index));
			}

			int forestDepth = Integer.parseInt(leftAndRightCountStringConverter
					.toString());

			System.out.println(leftCount);
			System.out.println(forestDepth);
			// END OF LEFT AND RIGHT INITIAL NODE PARSE BLOCK

			forestInfo = new ForestInfo(id, name, leftCount, forestDepth);
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
			while ((sCurrentLine = br.readLine()) != null) { // error checking
				System.out.println("Record:\t" + sCurrentLine);
				if (!sCurrentLine.startsWith("a"))
					break; // break as its the last line of the file.

				// ID PARSE BLOCK
				int index = 2;
				StringBuilder leafInfoStringBuilder = new StringBuilder();
				leafInfoStringBuilder.append(sCurrentLine.charAt(index));
				while (sCurrentLine.charAt(index + 1) != ' ') {
					index++;
					leafInfoStringBuilder.append(sCurrentLine.charAt(index));
				}
				int id = Integer.parseInt(leafInfoStringBuilder.toString());
				leafInfoStringBuilder.delete(0, leafInfoStringBuilder.length());
				index += 2;

				// LEVEL PARSE BLOCK
				leafInfoStringBuilder.append(sCurrentLine.charAt(index));
				while (sCurrentLine.charAt(index + 1) != ' ') {
					index++;
					leafInfoStringBuilder.append(sCurrentLine.charAt(index));
				}
				int level = Integer.parseInt(leafInfoStringBuilder.toString());
				leafInfoStringBuilder.delete(0, leafInfoStringBuilder.length());

				index += 2;

				// ANC PARSE BLOCK
				leafInfoStringBuilder.append(sCurrentLine.charAt(index));
				while (sCurrentLine.length() != index + 1
						&& sCurrentLine.charAt(index + 1) != ' ') {
					index++;
					leafInfoStringBuilder.append(sCurrentLine.charAt(index));
				}
				int anc = Integer.parseInt(leafInfoStringBuilder.toString());

				// FINAL DEBUG AND ADDITION OF ANC AND LEVEL TO THEIR INFO
				// QUEUES
				System.out.println("Id, Level, ANC " + id + level + anc);
				leafInfoStringBuilder.delete(0, leafInfoStringBuilder.length());

				info.addAnc(anc);
				info.addLevel(level);

			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		br.close();
		return info;
	}
}
