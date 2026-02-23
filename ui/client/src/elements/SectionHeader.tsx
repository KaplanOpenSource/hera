import { Typography } from "@mui/material";

export const SectionHeader = ({
  children,
}: {
  children: React.ReactNode;
}) => (
  <Typography variant="overline" color="text.secondary" sx={{ mt: 2, mb: 1, display: "block" }}>
    {children}
  </Typography>
);
