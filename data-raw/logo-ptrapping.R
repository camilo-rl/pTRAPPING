library(hexSticker)
library(showtext)

img <- "data-raw/ptrapping_logo.png"

## Loading Google fonts (http://www.google.com/fonts)

font_add_google("Comic Neue", "comic")

sticker(
  subplot = img,
  package = "pTRAPPING",
  p_color = "white",
  p_size = 17.8,
  p_y = 0.5,
  p_family = "comic",
  p_fontface = "bold",
  filename = "man/figures/logo.png",
  s_x = 0.99,
  s_y = 1,
  s_width = .8,
  h_fill = "#171717",
  h_color = "#000000"
)
