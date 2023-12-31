using System.IO;
using System.Reflection;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace VADER
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        private void HelloButton_Checked(object sender, RoutedEventArgs e)
        {

        }

        private void GoodbyeButton_Checked(object sender, RoutedEventArgs e)
        {

        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            if (HelloButton.IsChecked == true)
            {
                MessageBox.Show("Hello.");
            }
            else if (GoodbyeButton.IsChecked == true)
            {
                MessageBox.Show("Goodbye.");
            }
        }

        private void Button_Click_WriteToFile(object sender, RoutedEventArgs e)
        {
            
            string currentDirectory = Directory.GetCurrentDirectory();
            using (StreamWriter outputFile = new StreamWriter(System.IO.Path.Combine(currentDirectory, OutputFilenameTextBox.Text)))
            {
                string[] lines;
                if (HelloButton.IsChecked == true)
                {
                    lines = new string[] { "Hello" };
                }
                else if (GoodbyeButton.IsChecked == true)
                {
                    lines = new string[] { "Goodbye" };
                } else
                {
                    lines = new string[] { "Error" };
                }
                foreach (string line in lines)
                    outputFile.WriteLine(line);
            }
        }
    }
}