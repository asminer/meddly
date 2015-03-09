package logging;

import java.io.FileOutputStream;
import java.io.PrintStream;

public class Logger {
	PrintStream fileToWrite;

	public Logger(PrintStream fileToWrite) {
		this.fileToWrite = fileToWrite;
	}

	public Logger() {
		this.fileToWrite = System.out;
	}

	public void log(String Message) {
		fileToWrite.println(Message);
	}
}
