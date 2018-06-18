# ggdc-submitter

## Automatically downloading GGDC results with Microsoft Outlook

source: https://www.pixelchef.net/content/rule-autosave-attachment-outlook
VBA Code for Outlook:
```
Public Sub saveAttachtoDisk (itm As Outlook.MailItem)
Dim objAtt As Outlook.Attachment
Dim saveFolder As String
saveFolder = "c:\temp\"
  For Each objAtt In itm.Attachments        
    if InStr(objAtt.DisplayName, '.csv') Then
      objAtt.SaveAsFile saveFolder & "\" & objAtt.DisplayName
        end if
  Next
End Sub
```
