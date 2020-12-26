#include <QtWidgets/QApplication>
#include"ModelShader.h"
int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	ModelShader w(NULL);
	w.show();
	w.setModel("./3DModelFile/Moon2K.obj");

	return a.exec();
}
 