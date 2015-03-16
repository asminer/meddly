package main;

import info.LeafInfo;

import java.io.IOException;

import javafx.animation.Animation;
import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.application.Application;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.Scene;
import javafx.scene.chart.BarChart;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.stage.Stage;
import javafx.util.Duration;
import logic.ForestInfoParser;

public class ApplicationExecution extends Application {
	private LeafInfo info;
	private static ForestInfoParser applicationInfoParser = new ForestInfoParser();

	@SuppressWarnings("unchecked")
	@Override
	public void start(Stage stage) {
		final NumberAxis xAxis = new NumberAxis();
		final CategoryAxis yAxis = new CategoryAxis();
		final BarChart<Number, String> bc = new BarChart<Number, String>(xAxis,
				yAxis);
		Timeline tl = new Timeline();
		bc.setTitle("Summary of Forest Count");
		bc.setAnimated(true);
		xAxis.setLabel("Number of nodes");
		xAxis.setTickLabelRotation(90);
		yAxis.setLabel("Forest Level");
		XYChart.Series series = null;

		try {
			
			series = applicationInfoParser.initalizeForestInfo();
			info = applicationInfoParser.readNodeInfoFromJsonFile();
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (org.json.simple.parser.ParseException e) {
			e.printStackTrace();
		}

		
		tl.getKeyFrames().addAll(
				new KeyFrame(Duration.millis(200),
						new EventHandler<ActionEvent>() {
							@Override
							public void handle(ActionEvent actionEvent) {
								XYChart.Series<Number, String> series = bc
										.getData().get(0);
								try {
									if(info.hasNext()) {
										int level = info.getLevel();
										int anc = info.getAnc();
										
										XYChart.Data<Number, String> x = (XYChart.Data<Number, String>) series
												.getData()
												.get((int) (level - 1));
										System.out.print("Grabbed node at position: " + (level - 1) );
										System.out.println(" Value at node is: " + (x.getXValue().longValue()));
										System.out.print("Updating value anc: " + anc);
										x.setXValue(x.getXValue().intValue() + anc);
										System.out.println(" Value at node is: " + (x.getXValue().longValue()));
										
									}

								} catch (Exception e) {
									e.printStackTrace();
								}

							}
						}));
		tl.setCycleCount(Animation.INDEFINITE);
		tl.play();
		Scene scene = new Scene(bc, 800, 600);
		bc.getData().addAll(series);
		stage.setScene(scene);
		stage.show();
	}

	
	
	
	
	public static void main(String[] args)
			throws org.json.simple.parser.ParseException {

		launch(args);

	}
}
