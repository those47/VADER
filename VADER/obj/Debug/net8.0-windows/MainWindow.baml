   M S B A M L   `   `   `  ����@  <VADER, Version=1.0.0.0, Culture=neutral, PublicKeyToken=nullclr-namespace:VADERVADER  5   	       VADER.MainWindow   5   	   W SPresentationCore, Version=8.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35\ XPresentationFramework, Version=8.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35e aSystem.Windows.Controls.Ribbon, Version=8.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089R NWindowsBase, Version=8.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35F 9http://schemas.microsoft.com/winfx/2006/xaml/presentation     5   	   R NSystem.Xaml, Version=8.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089:x,http://schemas.microsoft.com/winfx/2006/xaml     5   	   @d2http://schemas.microsoft.com/expression/blend/2008     5   	   Jmc;http://schemas.openxmlformats.org/markup-compatibility/2006     5   	   %localclr-namespace:VADER     5   	     � Title$  
MainWindow��5   	   $	��450��6   $	��800��6)   .��5	      � .+�5
   
   �� $
��Left=�6   $��535,274,0,0q�6/   $
��Wrap=�6D   $��Output File:��6X   $	��Top=�6l   	� 5   
   -   5         	�Name$ HelloButton��5      $��Hello��6+    d� 	IsChecked$
 TrueY�6;   $
��Left=�6L   $��543,188,0,0q�6g   $	��Top=�6|   6�   	� 5   
   -   5        $ GoodbyeButton��5      $��	Goodbye
��6-   $
��Left=�6I   $��543,217,0,0q�6d   $	��Top=�6y   $��20��6�   6�   �� 5   
   -   5        $��Display��5      $��Center=�6#   $��	0,325,0,0q�6@   $	��Top=�6S   6k   �� 5   
   -   5         ��Name$ OutputFilenameTextBox��5      $
��Left=�61   $��506,306,0,0q�6L    �� TextWrapping$
 Wrap=�6a   $��
output.csv��6u   $	��Top=�6�   $	��120��6�   �� 5   
   $��Write to File��6   $��Center=�6)   $��	0,382,0,0q�6F   $	��Top=�6Y   � 5   
   $��45,33,479,239q�6    � ColumnDefinitions 5      �� 5      6%   �� 5      6%    � RowDefinitions 5      �� 5      6"   �� 5      6"   �� 5      6"   �� 5      6"   �� 5      6"   .+�5      �� $��0��6   $��0��6'   $��Center=�68   $
��Wrap=�6U   $��	Parameter��6i   $��Center=�6z   �� 5      $��0��6   $��1��6'   $��Center=�68   $
��Wrap=�6U   $��Input��6i   $��Center=�6v   �� 5      $��1��6   $��0��6'   $��Center=�68   $
��Wrap=�6U   $��Stream Length (m)��6i   $��Center=�6�   �� 5      $��1��6   $��1��6%   $��Center=�67   $
 Wrap=�6T   $��50��6h   $��Center=�6r   �� 5       $��2��6   $��0��6'   $��Center=�68   $
��Wrap=�6U   $��Segment Length (m)��6i   $��Center=�6�   �� 5!      $��2��6   $��1��6%   $��Center=�67   $
 Wrap=�6T   $	��0.1��6h   $��Center=�6s   �� 5"      $��3��6   $��Center=�6(   $
��Wrap=�6E   $��Duration (s)��6Y   $��Center=�6m   �� 5#      $��3��6   $��1��6%   $��Center=�67   $
 Wrap=�6T   $
��3600��6h   $��Center=�6t   �� 5$      $��4��6   $��0��6'   $��Center=�68   $
��Wrap=�6U   $��Simulation Interval (s)��6i   $��Center=�6�   �� 5%      $��4��6   $��1��6%   $��Center=�67   $
 Wrap=�6T   $��1��6h   $��Center=�6q   5&      5(      5)      