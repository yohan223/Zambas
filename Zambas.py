import sys
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QVBoxLayout, QLabel, QFileDialog, QTextEdit
from pyproj import Transformer

class ProjectionApp(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()
        
        self.label = QLabel("Projection Transformation Tool", self)
        layout.addWidget(self.label)
        
        self.input_text = QTextEdit(self)
        self.input_text.setPlaceholderText("Enter coordinates here...")
        layout.addWidget(self.input_text)
        
        self.transform_button = QPushButton("Transform Coordinates", self)
        self.transform_button.clicked.connect(self.transform_coordinates)
        layout.addWidget(self.transform_button)
        
        self.output_text = QTextEdit(self)
        self.output_text.setPlaceholderText("Output appears here...")
        layout.addWidget(self.output_text)
        
        self.setLayout(layout)
        self.setWindowTitle("Projection Transformation")

    def transform_coordinates(self):
        # Example transformation: WGS 84 to UTM Zone 33S
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:32733", always_xy=True)
        input_data = self.input_text.toPlainText().strip()
        
        if input_data:
            try:
                coords = [float(x) for x in input_data.split()]
                if len(coords) == 2:
                    x, y = transformer.transform(coords[0], coords[1])
                    self.output_text.setText(f"Transformed: {x}, {y}")
                else:
                    self.output_text.setText("Error: Enter two numeric values separated by space.")
            except ValueError:
                self.output_text.setText("Error: Invalid input format.")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ProjectionApp()
    window.show()
    sys.exit(app.exec_())

