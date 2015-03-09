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
		bc.setTitle("Summary");
		xAxis.setLabel("Number of nodes");
		xAxis.setTickLabelRotation(90);
		yAxis.setLabel("Forest Level");
		XYChart.Series series = null;

		try {
			series = applicationInfoParser.initalizeForestInfoFromJsonFile();
			info = applicationInfoParser.readNodeInfoFromJsonFile();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (org.json.simple.parser.ParseException e) {
			e.printStackTrace();
		}

		Timeline tl = new Timeline();
		tl.getKeyFrames().add(
				new KeyFrame(Duration.millis(100),
						new EventHandler<ActionEvent>() {
							@Override
							public void handle(ActionEvent actionEvent) {
								XYChart.Series<Number, String> series = bc
										.getData().get(0);
								try {
									while (info.hasNext()) {
										XYChart.Data<Number, String> x = (XYChart.Data<Number, String>) series
												.getData()
												.get((int) (info.getLevel() - 1));
										x.setXValue(x.getXValue().longValue()
												+ info.getAnc());
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
