// Give the browser time to start the download before the bytes are freed.
const REVOKE_DELAY_MS = 10000;

// Makes the browser save the bytes as a file on the user's machine.
export const saveFile = (bytes: ArrayBuffer, fileName: string) => {
  const url = URL.createObjectURL(new Blob([bytes], { type: 'application/zip' }));
  const link = document.createElement('a');
  link.href = url;
  link.download = fileName;
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), REVOKE_DELAY_MS);
}
