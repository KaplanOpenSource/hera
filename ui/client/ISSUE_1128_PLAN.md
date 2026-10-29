# Issue 1128 - Browse the project folder

Click the folder path under the project in the tree. A tab opens on the right
and shows the folder as a tree. You can go down into folders, never above the
project folder.

The folder is `desc.filesDirectory` in the project config. It can be edited,
empty, or missing from disk.

## Steps

1. Server (Python - needs your OK): a route that takes a project name and a
   path from the project folder. It returns nothing if the project has no
   folder, the path goes above the project folder, or the path does not exist.
   Else it returns the file names in it.
2. Client: clicking the path opens a "Files" tab. The copy button stays.
3. Client: the tab shows an MUI tree. It never reads the whole tree, which can
   be big. At first it shows only the project folder's contents. Clicking a
   subfolder asks the server for its files and subfolders, then shows them.
   It shows a message if the folder is not found.
4. Tests for the route, the click, and loading a folder.
